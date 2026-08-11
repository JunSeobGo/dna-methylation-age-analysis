# age-density 후보 모델 잠금과 새 외부 검증 계획

## 현재 결론

중첩 교차검증에서 선정한 `age_density`를 8개 학습 코호트 5,541명과 869개 CpG만으로 잠갔습니다. 모델 선택·결측 대치·가중치 계산에는 외부 검증 데이터를 사용하지 않았습니다.

이 모델은 아직 배포 모델이 아니라 **외부 검증 대기 후보**입니다. 기존 `GSE87571`은 표준 모델의 평가에 이미 사용했으므로 새 후보 평가에 재사용하지 않습니다.

## 잠금 결과

| 항목 | 값 |
| --- | ---: |
| 가중치 | 연령 구간 역빈도 기반 `age_density`, 개별 상한 5 |
| 학습 표본 | 5,541명, 8개 코호트 |
| 입력 특징 | 869 CpG |
| 학습 연령 범위 | 1~103세 |
| 선택 alpha | 0.05 |
| 선택 lambda | 0.7471671 |
| 비영 CpG | 377개 |
| 학습 전용 LOCO tuning MAE | 3.3776년 |
| 중첩 CV 후보 MAE | 3.7810년 |
| 외부 검증 사용 | 없음 |

LOCO tuning MAE는 같은 학습 코호트에서 하이퍼파라미터를 선택하기 위한 값이므로 일반화 성능으로 보고하지 않습니다. 후보의 내부 성능 주장은 하이퍼파라미터 선택과 평가가 분리된 중첩 CV MAE 3.7810년을 사용합니다.

잠금 아티팩트는 로컬에 다음과 같이 생성됩니다. 대용량·재생성 가능 파일이므로 Git에는 커밋하지 않습니다.

- `data/processed/locked_age_density_candidate.rds`: glmnet 객체, CpG 순서, 학습 중앙값, 가중치 정책
- `outputs/age_density_candidate_lock_manifest.csv`: 모델·학습 bundle MD5, 패키지 버전, 잠금 시각
- `outputs/age_density_candidate_selection.csv`: 선택 파라미터와 학습 범위
- `outputs/age_density_candidate_nonzero_coefficients.csv`: 비영 계수 감사표

## 새 외부 검증 후보

### 1. 주 검증: SATSA `S-BSST1206`

[EMBL-EBI BioStudies의 SATSA](https://www.ebi.ac.uk/biostudies/studies/S-BSST1206)는 혈액 450K·EPIC beta와 표현형을 공개합니다. 고령층 표본이 충분해 현재 가장 중요한 80세 이상 과소예측을 검증하는 주 코호트로 선정했습니다.

- 전체 표현형: 1,469행
- 450K: 1,094행, 447명
- 반복측정: 개인당 최대 5회
- 쌍둥이 가족: 266개 `TID` cluster
- 개인별 가장 이른 450K 표본 1개 선택: 447명, 48.027~97.879세
- 선택 표본 중 80세 이상 81명, 90세 이상 9명
- 성별: 여성 265명, 남성 182명
- 학습 표본 ID와 문자상 중복: 0건

반복측정 전체를 독립 표본처럼 사용하면 신뢰구간이 과도하게 좁아지므로, 주 분석은 개인별 최초 표본만 사용합니다. 쌍둥이 의존성은 `TID` 단위 paired cluster bootstrap으로 반영합니다.

표현형 manifest는 모델을 먼저 잠근 뒤 생성했습니다. beta 5개 파일의 공식 URL·바이트와 스트리밍 추출 품질 게이트도 모델 평가 전에 고정했습니다. 자세한 실행 방법은 [SATSA beta 스트리밍 추출과 품질 게이트](satsa-beta-extraction.md)를 참고합니다.

### 2. 보조 검증: `GSE193879` 건강 대조군

[GSE193879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193879)은 EPIC whole blood 127명 중 비감염 건강 대조군 69명을 포함합니다. 약 0.1~17세 범위를 제공하므로 20세 미만 편향의 보조 외부 검증 후보입니다.

MIS-C와 COVID-19 표본은 평가에서 제외하고, 건강 대조군의 sample-level 질환 라벨과 연령을 beta 다운로드 전에 먼저 잠가야 합니다. 질환 연구에서 모집된 대조군이므로 결과는 소아 일반집단 검증이 아니라 소아 민감도 분석으로 해석합니다.

### 3. 보조 검증: `GSE201287` 건강 대조군

[GSE201287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201287)은 450K whole blood 80명 중 건강 대조군 40명을 포함합니다. 중국인 표본에서 인구집단 이동을 점검할 수 있지만 표본이 작고 MDD 환자와 섞여 있으므로 건강 대조군만 분리한 보조 분석으로 제한합니다.

### 제외 후보

`GSE85210`과 `GSE67393`은 혈액 450K와 공개 beta 조건은 충족하지만 공개 GEO 표본 메타데이터에 연령이 없어 연령 예측 외부 검증에서 제외했습니다. beta가 있다는 이유만으로 연령 라벨이 없는 코호트를 내려받지 않습니다.

## 외부 검증 전에 고정한 판정 기준

SATSA에서 `age_density`와 기존 표준 잠금 모델을 동일 표본에 적용해 paired comparison을 수행합니다.

1. 전체 MAE가 표준 모델보다 0.20년 넘게 악화되지 않아야 합니다.
2. 80세 이상 MAE가 최소 0.50년 개선되어야 합니다.
3. 절대 평균 편향이 표준 모델보다 0.50년 넘게 악화되지 않아야 합니다.
4. `TID` cluster bootstrap에서 80세 이상 MAE 개선의 95% 신뢰구간 하한이 0보다 커야 합니다.

소아 보조 검증은 `GSE193879` 건강 대조군에서 전체 MAE 악화 0.20년 이내, 20세 미만 MAE 개선 0.50년 이상을 확인합니다. SATSA와 소아 검증을 모두 통과하기 전에는 기존 표준 모델을 교체하지 않습니다.

## 실행 방법

후보 잠금 설계만 확인합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/14_lock_age_density_candidate.R --dry-run
```

학습 데이터만으로 후보를 잠급니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/14_lock_age_density_candidate.R --confirm-lock
```

SATSA 표본 규칙을 확인하고 공개 표현형만 내려받아 manifest를 잠급니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/15_prepare_satsa_external_manifest.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/15_prepare_satsa_external_manifest.R --download-metadata
```

## 다음 실행 단계

1. ~~SATSA beta 5개 파일을 내려받아 공식 바이트와 로컬 MD5를 확인합니다.~~ 완료
2. ~~전체 beta를 메모리에 올리지 않고 869개 CpG와 잠긴 447개 표본만 스트리밍 추출합니다.~~ 완료
3. ~~사전 고정한 feature coverage, beta 범위, 결측률, 표본 순서·중복 품질 게이트를 판정합니다.~~ 완료, 판정 `FAIL`
4. 품질 게이트 통과 시에만 잠긴 두 모델을 한 번 적용하고 사전 판정 기준 및 `TID` cluster bootstrap을 계산합니다. **미수행**
5. 이후에만 `GSE193879` 건강 대조군의 소아 민감도 분석을 진행합니다. **미수행**

3단계에서 869개 모델 CpG 중 610개(70.20%)만 SATSA 공개 행렬에 존재해 사전 기준 0.95에 미달했습니다. 사전에 정한 대로 누락 CpG를 대치하거나 기준을 완화하지 않았고, 4단계 이후를 진행하지 않았습니다. 상세 수치는 [SATSA beta 스트리밍 추출과 품질 게이트](satsa-beta-extraction.md)를 참고합니다.

[ADR-0007](adr/0007-satsa-coverage-failure-response.md)에 따라 이 문서가 다루는 **869 CpG 트랙의 SATSA 주 검증은 수행 불가로 종결**했습니다. 이 문서의 잠금 결과와 판정 기준은 기록으로 보존하며 수정하지 않습니다. 재검증은 SATSA 가용 610개 CpG로 새로 잠근 별도 모델 쌍에서 진행합니다.
