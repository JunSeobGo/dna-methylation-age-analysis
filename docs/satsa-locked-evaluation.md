# SATSA 잠금 모델 비교 평가

## 목적

SATSA 외부 beta 결과를 확인하기 전에 표준 잠금 모델과 `age_density` 후보의 비교 방법을 실행 가능한 코드와 설정으로 고정합니다. 이 문서의 기준은 SATSA 성능을 확인한 뒤 변경하지 않습니다.

## 실행 전제

다음 입력이 모두 잠금 상태와 일치해야 합니다.

- 표준 잠금 모델: `data/processed/locked_elastic_net_model.rds`
- age-density 후보: `data/processed/locked_age_density_candidate.rds`
- 품질 게이트 통과 bundle: `data/processed/satsa_external_validation_bundle.rds`
- 표본: 개인별 가장 이른 450K 표본 447명
- 입력 특징: 두 모델에 공통인 동일 순서의 869 CpG
- 80세 이상: 사전 잠금 manifest와 같은 81명
- 후보·표준 모델 MD5: bundle에 기록된 값과 실제 파일이 일치

하나라도 다르면 성능을 계산하지 않습니다.

## 기계 판독 평가 정책

`config/satsa_evaluation_policy.csv`가 다음 기준과 bootstrap 설정을 고정합니다.

| 조건 | 통과 기준 |
| --- | ---: |
| 전체 MAE 비열등성 | 후보 MAE - 표준 MAE ≤ 0.20년 |
| 80세 이상 MAE 개선 | 표준 MAE - 후보 MAE ≥ 0.50년 |
| 전체 절대 평균 편향 악화 | 후보 절대 bias - 표준 절대 bias ≤ 0.50년 |
| 80세 이상 불확실성 | MAE 개선의 TID bootstrap 95% CI 하한 > 0 |

네 조건을 모두 통과할 때만 `replace_standard_model = TRUE`가 됩니다. 하나라도 실패하면 표준 모델을 유지합니다.

## 가족 단위 불확실성

SATSA는 쌍둥이 코호트이므로 80세 이상 표본을 개인 단위로 독립 재표집하지 않습니다.

1. 80세 이상 표본의 고유 `TID` 목록을 만듭니다.
2. 같은 수의 `TID`를 복원추출합니다.
3. 선택된 `TID`에 속한 쌍둥이 표본을 함께 포함합니다.
4. 동일 표본에서 `표준 MAE - 후보 MAE`를 계산합니다.
5. seed `20260804`로 10,000회 반복하고 95% 신뢰구간을 계산합니다.

두 모델이 항상 같은 bootstrap 표본에 적용되므로 paired comparison입니다.

## 실행 보호

외부 bundle을 읽지 않고 잠금 모델·정책만 확인합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/17_evaluate_satsa_locked_models.R --dry-run
```

SATSA 추출 품질 게이트가 통과한 뒤 최초 평가를 한 번 실행합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/17_evaluate_satsa_locked_models.R --confirm-external-evaluation
```

기존 결과 파일이 하나라도 있으면 기본적으로 재실행하지 않습니다. 입력 오류 등 정당한 사유가 있을 때만 이유를 함께 기록합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/17_evaluate_satsa_locked_models.R --confirm-external-evaluation --force "--rerun-reason=입력 파일 오류 수정 후 사전 정책 그대로 재평가"
```

## 로컬 산출물

- `outputs/satsa_locked_metrics.csv`: 모델별 전체·80세 이상 지표
- `outputs/satsa_locked_predictions.csv`: 동일 표본의 두 모델 예측과 paired 오차 개선
- `outputs/satsa_locked_tid_bootstrap.csv`: 10,000회 TID bootstrap 개선값
- `outputs/satsa_locked_decision.csv`: 네 조건과 최종 교체 판정
- `outputs/satsa_locked_evaluation_manifest.csv`: 정책·bundle·모델 MD5와 실행 시각
- `outputs/satsa_locked_diagnostics.png`: 두 모델의 예측·잔차 진단

모든 산출물은 재생성 가능하거나 외부 평가 결과를 포함하므로 Git에서 제외합니다. 최초 결과의 핵심 수치와 해석만 평가 완료 후 문서에 추가합니다.
