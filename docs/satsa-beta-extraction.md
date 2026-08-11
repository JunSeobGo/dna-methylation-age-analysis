# SATSA beta 스트리밍 추출과 품질 게이트

## 목적

SATSA `S-BSST1206`의 공개 beta 행렬에서 사전에 잠근 447명과 두 모델이 공통으로 사용하는 869개 CpG만 추출합니다. 이 단계는 입력 호환성을 검증하고 외부 검증 bundle을 만드는 단계이며, 모델 예측이나 성능 비교는 수행하지 않습니다.

## 공식 원본 구조

[EMBL-EBI BioStudies의 SATSA](https://www.ebi.ac.uk/biostudies/studies/S-BSST1206)는 `SATSA_DNA_meth_1.txt`부터 `SATSA_DNA_meth_5.txt`까지 5개 파일을 제공합니다.

- 전체 크기: 6,743,832,915바이트(약 6.74GB)
- 전체 행: 250,816 CpG
- 열: 1,469개 SATSA SID
- 구성: 450K 1,094개와 EPIC 375개를 함께 전처리한 beta
- 파일당 약 5만 행이며 모든 파일의 SID header 순서는 동일해야 함

URL과 공식 바이트 크기는 `config/dataset_download_manifest.csv`에 고정했습니다. 다운로드가 끝나면 로컬 MD5를 계산해 `outputs/satsa_beta_file_audit.csv`에 남깁니다.

## 메모리 절약 방식

`16_prepare_satsa_external_bundle.R`은 전체 6.74GB 행렬을 한 번에 읽지 않습니다.

1. 파일 header에서 잠근 447개 SID의 열 위치를 찾습니다.
2. 한 번에 1,000행씩 읽어 모델의 869개 CpG만 선택합니다.
3. 선택된 행에서 447개 beta만 숫자로 변환합니다.
4. 파일 하나가 끝날 때마다 로컬 checkpoint를 저장합니다.

최종 행렬 크기는 `447 × 869`이며, 원본·checkpoint·bundle·감사 산출물은 모두 Git에서 제외됩니다.

## 데이터 확인 전 사전 고정 품질 게이트

다음 기준은 SATSA 모델 성능을 보거나 beta 추출 결과를 확인하기 전에 고정했습니다.

- 5개 파일의 실제 바이트가 공식 manifest와 정확히 일치해야 함
- 5개 파일의 총 CpG 행 수가 250,816개여야 함
- 각 파일 header에 1,469개의 고유 SID가 있어야 함
- 잠근 447개 SID가 모든 파일 header에 정확히 존재해야 함
- 모델 CpG 869개 중 95% 이상이 SATSA에 존재해야 함
- 각 선택 표본의 결측률이 5% 이하여야 함
- 관측 beta 값이 모두 0~1 범위여야 함
- 후보 모델과 표준 모델의 869개 CpG 이름과 순서가 동일해야 함

게이트 하나라도 실패하면 감사 CSV는 남기되 외부 검증 bundle과 성능 결과는 만들지 않습니다. 누락 CpG를 결과에 맞춰 임의 대체하거나 품질 기준을 사후 완화하지 않습니다.

## 실행 방법

원본을 읽거나 내려받지 않고 잠금 상태와 파일 현황만 확인합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/16_prepare_satsa_external_bundle.R --dry-run
```

중단 시 `.part` 파일부터 이어받을 수 있도록 5개 원본을 다운로드합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/16_prepare_satsa_external_bundle.R --download
```

전송이 중단되면 같은 명령을 다시 실행합니다. `--dry-run`은 완성 파일뿐 아니라 `.part`의 현재 바이트와 진행률도 표시합니다. `.part`를 수동으로 완성 파일 이름으로 바꾸면 공식 크기 검사를 우회할 수 있으므로 이름은 직접 바꾸지 않습니다.

모델 성능을 계산하지 않고 스트리밍 추출과 품질검사만 수행합니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/16_prepare_satsa_external_bundle.R --extract --confirm-extraction
```

다운로드와 추출을 연속으로 실행할 수도 있습니다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/16_prepare_satsa_external_bundle.R --download --extract --confirm-extraction
```

## 로컬 산출물

- `data/processed/satsa_beta_extraction_checkpoint.rds`: 파일별 재개 checkpoint
- `data/processed/satsa_external_validation_bundle.rds`: 품질 게이트 통과 시에만 생성
- `outputs/satsa_beta_file_audit.csv`: 크기, MD5, 행·열, SID 연결 감사
- `outputs/satsa_feature_coverage.csv`: 869개 CpG별 존재 여부와 결측률
- `outputs/satsa_beta_extraction_summary.csv`: 전체 품질 게이트 결과

## 실제 실행 결과: 품질 게이트 `FAIL`

5개 원본을 모두 내려받아 스트리밍 추출을 1회 실행했습니다. 파일 무결성과 표본 연결은 모두 통과했지만 **CpG coverage 기준에서 실패**했으므로 외부 검증 bundle을 만들지 않았고 모델 성능도 계산하지 않았습니다.

| 검사 | 기준 | 실제 | 판정 |
| --- | --- | ---: | --- |
| 5개 파일 바이트 일치 | 공식 manifest와 정확히 일치 | 6,743,832,915바이트 | 통과 |
| 총 CpG 행 수 | 250,816 | 250,816 | 통과 |
| header SID 수 | 1,469 | 1,469 | 통과 |
| 잠근 447 SID 연결 | 5개 파일 모두 447 | 447 | 통과 |
| 모델 CpG coverage | ≥ 0.95 | **0.7020 (869개 중 610개)** | **실패** |
| 표본별 결측률 | ≤ 0.05 | **0.2980 (모든 표본 동일)** | **실패** |
| 관측 beta 범위 | 0~1 | 0.0000813~0.9976 | 통과 |
| 선택 SID 중복 | 0 | 0 | 통과 |

파일별 추출 CpG는 224·145·96·69·76개로 합계 610개이며 파일 간 중복은 없었습니다. 표본별 결측률 0.2980은 정확히 `259 / 869`로, 개별 표본의 측정 실패가 아니라 **없는 CpG 259개가 모든 표본에 공통으로 비어 있기 때문**입니다.

### 원인

SATSA 공개 행렬은 450K 전체 probe가 아니라 자체 품질검사 뒤 남은 250,816개 probe만 제공합니다. 학습 코호트 8개의 교집합으로 정한 869개 CpG 중 259개가 이 축소 집합에 포함되지 않았습니다. 코드나 표본 연결 문제가 아니라 **원본 probe 구성 차이**입니다.

### 영향

없는 259개는 후보 모델의 비영 계수에도 걸쳐 있습니다.

| 항목 | 값 |
| --- | ---: |
| `age_density` 후보의 비영 CpG | 377개 |
| 그중 SATSA에 없는 CpG | 108개 |
| 비영 계수 절댓값 합 | 849.74 |
| 없는 CpG가 차지하는 절댓값 합 | 224.77 (26.4%) |

잠금 모델의 학습 중앙값으로 259개를 대치하면 예측이 학습 평균 쪽으로 강하게 수축합니다. 이는 고령층 과소예측을 검증하려는 이번 평가의 목적 자체를 무효화하므로, 사전에 정한 대로 **게이트를 사후 완화하거나 누락 CpG를 대체하지 않았습니다.**

## 다음 단계

품질 게이트가 `PASS`일 때만 별도 스크립트에서 잠근 표준 모델과 `age_density` 후보를 한 번 적용합니다. 비교 기준과 `TID` cluster bootstrap 규칙은 기존 사전 계획을 그대로 사용합니다.

이번 실행은 `FAIL`이므로 `17_evaluate_satsa_locked_models.R`은 실행하지 않았고, 현재 869 CpG 잠금 모델 쌍으로는 SATSA 주 검증을 수행할 수 없습니다. 표준 모델을 교체하지 않는다는 기존 적용 경계는 그대로 유지됩니다.

후속 방향은 [ADR-0007](adr/0007-satsa-coverage-failure-response.md)에서 결정했습니다.

- **869 CpG 트랙**: SATSA 주 검증을 수행 불가로 종결합니다. 같은 데이터에 다른 기준을 적용해 재시도하지 않습니다.
- **610 CpG 제한 트랙**: SATSA에 실제로 존재하는 610개 CpG만 사용하는 별도 모델 쌍을 학습 데이터만으로 새로 잠그고, ADR-0006의 판정 기준을 값 변경 없이 그대로 적용해 1회 재검증합니다.

특징 목록은 SATSA beta의 존재 여부만으로 결정하며 연령 라벨과 예측 오차는 사용하지 않습니다. 두 트랙의 내부 성능은 특징 집합이 달라 직접 비교하지 않습니다.
