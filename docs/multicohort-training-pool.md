# 다중 코호트 학습 풀과 과적합 방지 정책

## 결론

표본 수를 늘리는 방향은 타당하지만 `GSE207605`의 11,910개 표본을 한꺼번에 무작위 분할하면 성능이 과대평가될 수 있다. 코호트마다 연령 분포, 질환 상태, 전처리 방식, 가족 구조가 다르기 때문이다.

따라서 19개 구성 코호트를 다음 세 정책으로 나눈다.

| 정책 | 코호트 수 | 표본 수 | 사용 방법 |
| --- | ---: | ---: | --- |
| `include` | 8 | 5,541 | 1차 다중 코호트 학습 후보 |
| `conditional` | 3 | 1,853 | 반복측정·암 상태·트라우마 영향을 통제할 수 있을 때만 민감도 분석에 추가 |
| `exclude` | 8 | 4,516 | 단일 연령, 질환 혼합, 기술·가족 의존성 때문에 1차 학습에서 제외 |

세부 근거와 원본 URL은 `config/gse207605_cohort_manifest.csv`에서 관리한다. 이 정책은 표본 수가 많은 코호트가 모델을 지배하는 것을 막기 위한 첫 번째 안전장치이며, 최종 확정은 원 메타데이터의 개체 식별자와 질환 라벨을 추가로 확인한 뒤 수행한다.

## 독립 외부 검증 세트

`GSE87571`은 `GSE207605`의 19개 구성 코호트에 포함되지 않는다. 이 코호트는 특징 선택, 결측 대치, 스케일링, 하이퍼파라미터 선택에 사용하지 않고 마지막 한 번의 외부 평가에만 사용한다.

`GSE84727`은 `GSE207605`에 포함되어 있고 조현병 환자와 대조군이 섞여 있으므로 이 프로젝트의 독립 외부 검증 세트로 사용하지 않는다.

## 과적합 방지 원칙

1. 행 단위 무작위 분할 대신 `source_series_id`를 그룹으로 사용하는 leave-one-cohort-out 평가를 우선한다.
2. 같은 사람의 반복측정, 쌍둥이, 기술 반복 표본은 반드시 같은 fold에 둔다.
3. CpG 선택, 결측 대치, 스케일링, 모델 튜닝은 각 학습 fold 안에서만 수행한다.
4. 코호트별 MAE, RMSE, 중앙절대오차와 연령대별 오차를 함께 보고한다.
5. `GSE87571` 외부 검증 결과를 확인하기 전까지 모델 선택을 끝낸다.
6. `GSE207605`의 CpG는 원 논문에서 연령 연관성으로 사전 선별된 집합이므로 탐색·벤치마크용으로 명시한다. 완전히 독립적인 특징 발견을 주장하지 않는다.

## 현재 확인된 구조

- 19개 파일의 총 표본 수: 11,910
- 코호트 간 중복 GSM ID: 0
- 1차 학습 후보 8개 코호트의 표본 수: 5,541
- 원본 파일별 CpG 수: 934~2,374개
- 1차 학습 후보에서 공통으로 관측되는 CpG만 모델 입력으로 사용

공통 CpG 수는 감사 스크립트가 실행 시점에 계산해 `data/processed/gse207605_primary_common_probes.txt`에 저장한다.

## 실행 방법

명세 형식만 점검하려면 다음을 실행한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/6_audit_gse207605_training_pool.R --dry-run
```

누락된 소형 재분석 파일을 내려받고 감사를 실행하려면 다음을 사용한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/6_audit_gse207605_training_pool.R --download
```

원본과 생성 결과는 Git에 포함하지 않는다. 정책 명세, 스크립트, 방법 문서만 버전 관리한다.

## 학습 bundle 생성

감사에서 확정한 1차 학습 정책을 바탕으로 `7_build_gse207605_training_bundle.R`이
Elastic Net 중첩 교차검증에 바로 투입할 수 있는 학습용 bundle을 만든다.

- `model_policy == include`인 8개 코호트만 선택한다.
- `data/processed/gse207605_primary_common_probes.txt`의 869개 공통 CpG만 동일 순서로 추출한다.
- 표본 방향(행=표본, 열=CpG)으로 정렬해 `5,541 × 869` beta 행렬을 만든다.
- `x`, `y`(chronological_age), `group`(source_series_id), `sample_id`, `feature_names`를
  담은 리스트를 `data/processed/gse207605_training_bundle.rds`로 저장한다.

결측 대치와 표준화는 교차검증의 각 학습 fold 안에서만 수행해야 하므로 이 단계에서는
전체 데이터 기준 대치를 하지 않고 `NA`를 그대로 보존한다. 생성 요약은
`outputs/gse207605_training_bundle_summary.csv`에 저장한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/7_build_gse207605_training_bundle.R
```

실행 시 다음을 검증한다.

- `X` 행 5,541 / 열 869, `y`·`group`·`sample_id` 길이 5,541
- `sample_id`·CpG 중복 0, beta 범위 0~1
- 코호트별 표본 수가 명세와 일치, 모든 행과 메타데이터 순서 일치
- 결측치 수와 비율 보고(대치하지 않음)

## 코호트 그룹 기반 중첩 교차검증 (Elastic Net 벤치마크)

`8_train_grouped_elastic_net.R`은 학습 bundle로 leave-one-cohort-out 중첩 교차검증을
수행한다. glmnet 패키지가 필요하다.

- **outer loop**: 8개 코호트를 하나씩 검증 세트로 두는 leave-one-cohort-out. outer test
  코호트는 inner tuning에 전혀 쓰지 않는다.
- **inner loop**: outer 학습 코호트 안에서 다시 leave-one-cohort-out 그룹 CV로 `alpha`,
  `lambda`를 선택한다.
- **누수 방지**: 결측 대치(열 중앙값)와 표준화는 각 학습 fold 안에서만 계산하고 검증·테스트
  fold에는 학습 fold 값만 적용한다. 전체 데이터 기준 대치·표준화를 하지 않는다.
- **가중치 비교**: `GSE55763`(2,639명)이 표본 수로 모델을 지배하지 않도록 일반 표본 가중치와
  코호트별 총 가중치를 같게 한 balanced 가중치를 함께 평가한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/8_train_grouped_elastic_net.R
```

### 벤치마크 결과 요약

전체 OOF 지표는 다음과 같다(수치는 재실행 시 데이터·환경에 따라 소폭 달라질 수 있다).

| 가중치 | MAE | RMSE | MedAE | R² |
| --- | ---: | ---: | ---: | ---: |
| standard | 약 4.15 | 약 5.18 | 약 3.57 | 약 0.91 |
| balanced | 약 4.47 | 약 5.52 | 약 3.94 | 약 0.90 |

- `standard`는 전체 MAE가 낮지만, 표본이 적고 연령이 극단인 코호트(초고령 `GSE30870`,
  소아 `GSE36054`)와 연령 구간 양끝(20세 미만, 80세 이상)에서 편향이 크다.
- `balanced`는 전체 MAE를 약간 희생하는 대신 소아·저연령·고연령 구간의 오차와 편향을 줄인다.
  이는 코호트 균형 가중치가 다수 코호트 지배를 완화한다는 것을 보여준다.
- 코호트별·연령 구간별 지표, fold별 `alpha`/`lambda`와 비영 계수 수, 예측-실제·잔차 진단
  그래프는 `outputs/gse207605_cv_*`에 저장된다(Git 제외).

> 869개 공통 CpG는 `GSE207605` 원 논문에서 연령 연관성으로 사전 선택된 특징이므로 이 결과는
> **벤치마크**로 해석한다. 최종 일반화 성능은 `GSE87571` 외부 검증에서 판단한다.

## GSE87571 외부 검증 bundle 준비

`9_prepare_gse87571_external_validation.R`은 최종 잠금 외부 검증 세트 `GSE87571`의 beta
행렬에서 869개 공통 CpG를 학습 bundle과 동일한 feature 순서로 추출해 외부 검증용 RDS를
만든다. **이 단계에서는 모델 성능(나이 기반 지표)을 계산하지 않는다.**

- `GSE87571_matrix1of2.txt.gz`, `GSE87571_matrix2of2.txt.gz`는 SWAN 정규화 Average Beta로,
  표본이 두 파일에 나뉘어 있고 각 표본은 두 열(`Xn`=beta, `Xn.1`=detection p-value)로 저장된다.
- beta 열은 exact-zero 비율로 자동 판정한다(detection p-value 열은 정확히 0인 값이 대부분).
- 열 라벨 `Xn`은 SOFT 메타데이터의 sample_title과 연결되고, 여기서 GSM과 age를 얻는다.
- 학습 데이터와 GSM이 겹치지 않는지, feature 순서가 학습 bundle과 정확히 일치하는지 검증한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/9_prepare_gse87571_external_validation.R
```

검증 결과 요약:

- 표본 732개 × CpG 869개, 연령 확인 가능 표본 729개(연령 14~94세)
- 학습 데이터와 GSM 중복 0건, feature 순서 학습 bundle과 일치
- 결측치는 대치하지 않고 보존, beta 범위 0~1
- 산출물: `data/processed/gse87571_external_validation_bundle.rds`,
  `outputs/gse87571_external_validation_summary.csv`(Git 제외)

## 다음 모델링 단계

1. 중첩 교차검증 결과를 바탕으로 코호트·CpG·대치 방법·`alpha`·`lambda`·가중치를 고정한다.
2. `conditional` 코호트를 하나씩 추가해 성능과 편향 변화를 민감도 분석한다.
3. 모델과 임계값을 고정한 뒤 `GSE87571`을 한 번만 평가한다.
