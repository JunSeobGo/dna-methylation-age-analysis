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
