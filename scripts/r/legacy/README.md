# 레거시 R 스크립트

프로젝트 초기 탐색 단계에서 사용한 기록입니다. **현행 모델 성능 재현에는 사용하지 않습니다.** 현행 파이프라인은 `scripts/r/`의 번호 순서 스크립트입니다.

| 파일 | 용도 |
| --- | --- |
| `Methyl_age_230526.R` | Horvath 21k CpG predictor로 DNAm age를 산출하던 초기 워크플로 |
| `Methyl_시각화(Ridge, Elastic net, Lasso).R` | 정규화 회귀 3종 비교와 시각화 실험 |
| `1_data_download_from_geo.r` | `1_data_download.R`과 중복인 GEO 다운로드 스크립트 |
| `GSM샘플 저장.R` | 표본 메타데이터 저장 보조 스크립트 |

## 실행 시 주의

이 스크립트들에는 과거 로컬 절대경로가 남아 있어 그대로 실행되지 않습니다. `setwd(...)`, `source(...)`, `read.csv(...)`, `read.csv.sql(...)` 경로를 현재 환경에 맞게 수정해야 합니다.

```r
project_dir <- normalizePath(".")
data_dir <- file.path(project_dir, "data", "raw")

source(file.path(data_dir, "AdditionalFile24NROMALIZATION.R.txt"))
dat0 <- sqldf::read.csv.sql(file.path(data_dir, "AdditionalFile26MethylationDataExample55.csv"))
```

`Methyl_age_230526.R`이 참조하는 입력 파일은 다음과 같습니다.

| 파일 | 역할 |
| --- | --- |
| `AdditionalFile21datMiniAnnotation27k.csv` | 27k CpG probe annotation |
| `AdditionalFile22probeAnnotation21kdatMethUsed.csv` | clock 계산에 필요한 CpG probe 목록 |
| `AdditionalFile23predictor.csv` | DNAm age predictor 계수 |
| `AdditionalFile24NROMALIZATION.R.txt` | 정규화 함수 |
| `AdditionalFile25StepwiseAnalysis.txt` | DNAm age 계산 로직 |
| `AdditionalFile26MethylationDataExample55.csv` | 예제 beta value 데이터 |
| `AdditionalFile27SampleAnnotationExample55.csv` | 예제 표본의 실제 연령 annotation |

입력 beta value CSV는 첫 번째 열에 Illumina CpG probe ID(예: `cg00000292`)를 두고 이후 열마다 표본의 numeric beta value를 포함해야 합니다. 필요한 probe가 없거나 숫자 형식이 아니면 계산이 중단되거나 `LogFile.txt`에 경고가 기록됩니다.

이 스크립트들은 `GEOquery`, `WGCNA`, `sqldf`, `readr`, `progress`에 의존하지만, 이 패키지들은 현행 파이프라인 재현의 필수 의존성이 아닙니다. 버전 정책은 `docs/r-environment-reproducibility.md`를 참고하세요.
