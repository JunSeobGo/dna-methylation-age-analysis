# R Scripts

This folder contains the R analysis scripts used for the project.

## 다중 코호트 학습 풀 감사

`6_audit_gse207605_training_pool.R`은 `GSE207605`의 19개 구성 코호트를 내려받고 표본 수, 연령 범위, beta 값 범위, GSM 중복, 1차 학습 후보의 공통 CpG를 검사한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/6_audit_gse207605_training_pool.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/6_audit_gse207605_training_pool.R --download
```

코호트별 포함·조건부·제외 판단은 `config/gse207605_cohort_manifest.csv`에서 관리하며, 상세 원칙은 `docs/multicohort-training-pool.md`를 참고한다.

## 1차 학습 bundle 생성

`7_build_gse207605_training_bundle.R`은 `model_policy == include`인 8개 코호트에서 869개 공통 CpG만 동일 순서로 추출해 `5,541 × 869` beta 행렬과 메타데이터를 `data/processed/gse207605_training_bundle.rds`로 저장한다. 결측치는 대치하지 않고 그대로 보존하며, 생성 요약은 `outputs/gse207605_training_bundle_summary.csv`에 남긴다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/7_build_gse207605_training_bundle.R
```

## 코호트 그룹 기반 중첩 교차검증

`8_train_grouped_elastic_net.R`은 학습 bundle로 leave-one-cohort-out 중첩 교차검증을 수행해 Elastic Net 벤치마크를 학습·평가한다(`glmnet` 필요). outer/inner loop를 모두 코호트 그룹 단위로 분리하고, 결측 대치와 표준화를 각 학습 fold 안에서만 계산해 데이터 누수를 막는다. 일반 가중치와 코호트 균형 가중치를 비교하며, 지표·fold 파라미터·진단 그래프는 `outputs/gse207605_cv_*`에 저장한다(Git 제외).

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/8_train_grouped_elastic_net.R
```

## GSE87571 외부 검증 bundle 준비

`9_prepare_gse87571_external_validation.R`은 최종 잠금 외부 검증 세트 `GSE87571`의 beta 행렬(`GSE87571_matrix1of2.txt.gz`, `GSE87571_matrix2of2.txt.gz`)에서 869개 공통 CpG를 학습 bundle과 동일한 feature 순서로 추출해 `data/processed/gse87571_external_validation_bundle.rds`로 저장한다. 표본은 두 파일에 나뉘어 있고 각 표본이 `Xn`(beta)와 `Xn.1`(detection p-value) 두 열로 저장되므로 beta 열을 exact-zero 비율로 판정한다. 학습 GSM과 중복이 없는지, feature 순서가 학습과 일치하는지 검증하며 **이 단계에서는 모델 성능을 계산하지 않는다.**

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/9_prepare_gse87571_external_validation.R
```

## 최종 모델 고정과 외부 검증

`10_lock_model_and_evaluate_gse87571.R`은 학습 코호트만으로 alpha, lambda, 결측 대치값, 가중치 방식을 고정한 뒤에만 `GSE87571` 외부 검증을 수행한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/10_lock_model_and_evaluate_gse87571.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/10_lock_model_and_evaluate_gse87571.R --confirm-external-evaluation
```

외부 평가 결과가 이미 있으면 재실행을 막으며, 명시적 재평가가 필요한 경우에만 `--force`를 사용한다.

## 학습곡선 측정

`11_measure_training_learning_curve.R`은 외부 검증셋을 제외하고 8개 학습 코호트만으로 표본 수별 성능을 측정한다. 코호트 하나를 통째로 홀드아웃한 상태에서 남은 코호트의 25%, 50%, 75%, 100%를 코호트·연령대별로 층화 추출한다. 고정 모델의 하이퍼파라미터를 그대로 사용해 표본 수 효과만 비교한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/11_measure_training_learning_curve.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/11_measure_training_learning_curve.R
```

설계, 결과 파일, 판정 기준은 `docs/learning-curve-analysis.md`에 정리했다.

## 연령·코호트 편향 감사

`12_audit_age_and_cohort_bias.R`은 학습곡선 100% 조건의 코호트 홀드아웃 예측과 학습 bundle을 사용해 연령 극단부 편향과 취약 코호트의 원인을 진단한다. sample_id·beta 범위·결측률을 확인하고, 각 홀드아웃 코호트에 대해 나머지 학습 코호트의 연령 지원도, CpG 중앙값 이동, 코호트별 MAE 부트스트랩 신뢰구간을 계산한다. `GSE87571`은 읽지 않는다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/12_audit_age_and_cohort_bias.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/12_audit_age_and_cohort_bias.R
```

해석과 다음 실험의 사전 채택 기준은 `docs/age-cohort-bias-audit.md`에 정리했다.

## 연령 편향 완화 후보 비교

`13_compare_age_bias_mitigation.R`은 외부 검증셋을 읽지 않고 8개 학습 코호트의 중첩 교차검증 안에서 표준·연령 밀도·연령/코호트 혼합 가중치와 inner OOF 선형 보정을 비교한다. 각 outer fold가 끝날 때 체크포인트를 저장하므로 실행이 중단돼도 완료 지점부터 다시 시작한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/13_compare_age_bias_mitigation.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/13_compare_age_bias_mitigation.R
```

결과와 채택 판단은 `docs/age-bias-mitigation-experiment.md`에 정리했다. 생성되는 `outputs/age_bias_mitigation_*` 파일은 재현 가능한 산출물이므로 Git에서 제외한다.

원본 beta 행렬은 다음에서 내려받는다(각 약 1.4GB, Git 제외).

```powershell
curl --fail --location --continue-at - --retry 8 --output data/raw/GSE87571/GSE87571_matrix1of2.txt.gz `
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/GSE87571_matrix1of2.txt.gz
curl --fail --location --continue-at - --retry 8 --output data/raw/GSE87571/GSE87571_matrix2of2.txt.gz `
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87571/suppl/GSE87571_matrix2of2.txt.gz
```

## Files
- 1_data_download.R: download GEO sample metadata
- 1_data_download_from_geo.r: duplicate download script retained for convenience
- GSM샘플 저장.R: sample metadata saving helper
- Methyl_age_230526.R: DNA methylation age modeling workflow
- Methyl_시각화(Ridge, Elastic net, Lasso).R: visualization and modeling experiments
