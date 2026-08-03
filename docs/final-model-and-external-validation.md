# 최종 모델 고정과 GSE87571 외부 검증

## 원칙

최종 Elastic Net의 가중치 방식, alpha, lambda, CpG별 결측 대치값은 `GSE207605` 학습 코호트 8개(5,541명)만으로 고정한다. 고정된 모델은 그 뒤에만 잠금 외부 검증 세트 `GSE87571`의 연령 라벨 표본 729개에 적용한다.

`GSE87571`은 모델 선택, CpG 선택, 결측 대치, 스케일링에 사용하지 않는다. 외부 beta 결측값도 학습 데이터에서 고정한 중앙값으로만 대치한다.

## 실행

외부 데이터를 사용하지 않고 선택 결과만 확인하려면 다음을 실행한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/10_lock_model_and_evaluate_gse87571.R --dry-run
```

최종 모델을 고정하고 외부 검증을 실행하려면 다음을 사용한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/10_lock_model_and_evaluate_gse87571.R --confirm-external-evaluation
```

결과 파일이 이미 있으면 잠금 외부 검증의 반복 실행을 막는다. 재실행이 불가피한 근거가 있을 때만 `--force`를 명시한다.

## 산출물

- `data/processed/locked_elastic_net_model.rds`: 모델, CpG 순서, 학습 중앙값, 고정 파라미터
- `outputs/final_model_selection.csv`: 학습 데이터만으로 선택한 설정
- `outputs/gse87571_external_metrics_overall.csv`: MAE, RMSE, MedAE, R², bias, calibration
- `outputs/gse87571_external_metrics_by_agebin.csv`: 연령 구간별 오차
- `outputs/gse87571_external_predictions.csv`: 표본별 예측과 잔차
- `outputs/gse87571_external_diagnostics.png`: 외부 진단 그림

원본 데이터와 생성 산출물은 Git에 포함하지 않는다.

## 최초 잠금 외부 검증 결과

2026-08-03에 수행한 최초 외부 검증에서는 내부 benchmark MAE가 더 낮았던 `standard` 가중치를 선택했다. 학습 코호트 leave-one-cohort-out 선택 결과는 alpha `0.05`, lambda `1.2626`, 비영 계수 295개였다.

| 평가 세트 | 표본 수 | MAE | RMSE | MedAE | R² |
| --- | ---: | ---: | ---: | ---: | ---: |
| GSE87571 | 729 | 2.80년 | 3.62년 | 2.34년 | 0.970 |

연령 구간별로는 80세 이상에서 MAE 6.08년, 평균 잔차 -5.74년으로 과소예측이 가장 컸다. 따라서 전체 MAE만으로 모델이 모든 연령대에 동일하게 정확하다고 해석하면 안 된다.

## 해석 경계

869개 CpG는 공개 재분석 자료에서 연령 연관성으로 사전 선별된 CpG 집합의 공통 부분이다. 따라서 결과는 이 특징 집합을 사용한 혈액 기반 연령 예측의 외부 검증으로 해석하며, 완전히 독립적인 CpG 발견 성능으로 주장하지 않는다.
