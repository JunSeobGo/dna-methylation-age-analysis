# DNA Methylation Age Analysis

DNA 메틸화 beta value를 이용해 DNA methylation age(DNAm age)를 산출하고, 실제 연령과의 차이(age acceleration)를 탐색하는 분석 프로젝트입니다. R 기반 DNAm age 산출 파이프라인과 Python/Jupyter 기반 전처리·탐색 분석을 함께 제공합니다.

> **연구용 분석 도구입니다.** 본 프로젝트의 결과는 임상적 진단, 치료 결정 또는 개인 건강 상태의 확정적 판단에 사용하면 안 됩니다.

## 분석 범위

- GEO(Gene Expression Omnibus)에서 GSM 샘플 메타데이터 수집
- Illumina CpG probe와 DNA 메틸화 beta value의 형식 검증
- 21k CpG probe 기반 DNAm age 계산
- 실제 연령과 예측 연령의 차이(`age_diff`) 분석
- BDI(Beck Depression Inventory) 점수와 age difference의 관계 탐색
- Ridge, Elastic Net, Lasso 회귀 실험 및 시각화

## 프로젝트 구조

```text
DNAage_clean/
├── data/
│   ├── raw/                 # 원천 데이터와 참조 파일
│   └── processed/           # 전처리 결과 저장 위치
├── docs/                    # 분석 배경 및 참고 문서
├── notebooks/               # Python 기반 전처리·탐색 분석 노트북
├── outputs/                 # 재생성 가능한 분석 결과
└── scripts/
    └── r/                   # R 기반 다운로드·DNAm age·모델링 스크립트
```

## 주요 구성 요소

| 위치 | 용도 |
| --- | --- |
| `scripts/r/1_data_download.R` | GEOquery로 지정한 GSM 샘플 메타데이터 다운로드 |
| `scripts/r/Methyl_age_230526.R` | CpG beta value에서 DNAm age 산출 및 결과 파일 생성 |
| `scripts/r/Methyl_시각화(Ridge, Elastic net, Lasso).R` | 정규화 회귀 모델 비교와 시각화 |
| `notebooks/2_preprocessing.ipynb` | 샘플 메타데이터 전처리 |
| `notebooks/4_BDI_and_AgeDiff_Correlation_Analysis.ipynb` | BDI와 `age_diff`의 상관·집단 비교 분석 |
| `notebooks/Methyl(파이썬 변환).ipynb` | DNAm age 파이프라인의 Python 변환 실험 |
| `data/raw/AdditionalFile2*.csv` | probe annotation, predictor, 예제 beta value 및 샘플 annotation |

## 요구 환경

### R

- R 4.2 이상 권장
- `BiocManager`, `GEOquery`, `WGCNA`, `sqldf`, `readr`, `progress`
- 다중 코호트 중첩 교차검증(`scripts/r/8_train_grouped_elastic_net.R`)에는 `glmnet` 필요

필요한 패키지는 아래와 같이 설치할 수 있습니다.

```r
install.packages(c("BiocManager", "WGCNA", "sqldf", "readr", "progress", "glmnet"))
BiocManager::install(c("GEOquery", "impute"))
```

### Python

- Python 3.10 이상 권장
- JupyterLab 또는 Jupyter Notebook
- `pandas`, `numpy`, `scikit-learn`, `scipy`, `matplotlib`, `seaborn`

```bash
python -m pip install pandas numpy scikit-learn scipy matplotlib seaborn jupyterlab
```

## 빠른 시작

### 1. 데이터 배치

분석에 필요한 참조 파일과 beta value 파일을 `data/raw/`에 둡니다. DNAm age 산출 스크립트가 사용하는 기본 입력 파일은 다음과 같습니다.

| 파일 | 역할 |
| --- | --- |
| `AdditionalFile21datMiniAnnotation27k.csv` | 27k CpG probe annotation |
| `AdditionalFile22probeAnnotation21kdatMethUsed.csv` | clock 계산에 필요한 CpG probe 목록 |
| `AdditionalFile23predictor.csv` | DNAm age predictor 계수 |
| `AdditionalFile24NROMALIZATION.R.txt` | 정규화 함수 |
| `AdditionalFile25StepwiseAnalysis.txt` | DNAm age 계산 로직 |
| `AdditionalFile26MethylationDataExample55.csv` | 예제 beta value 데이터 |
| `AdditionalFile27SampleAnnotationExample55.csv` | 예제 샘플의 실제 연령 annotation |

입력 beta value CSV는 **첫 번째 열에 Illumina CpG probe ID**(예: `cg00000292`)를 두고, 이후 열마다 샘플의 numeric beta value를 포함해야 합니다. 필요한 probe가 누락되었거나 숫자 형식이 아닌 값이 있으면 계산이 중단되거나 `LogFile.txt`에 경고가 기록됩니다.

### 2. R 스크립트의 경로 설정

`scripts/r/Methyl_age_230526.R`에는 과거 로컬 절대경로가 포함되어 있습니다. 실행 전 다음 항목을 현재 프로젝트 경로에 맞게 수정하세요.

- `setwd(...)`
- `source(...)`로 불러오는 `AdditionalFile24NROMALIZATION.R.txt`
- `read.csv(...)`, `read.csv.sql(...)`의 입력 경로

예시는 다음과 같습니다.

```r
project_dir <- normalizePath(".")
data_dir <- file.path(project_dir, "data", "raw")

source(file.path(data_dir, "AdditionalFile24NROMALIZATION.R.txt"))
dat0 <- sqldf::read.csv.sql(
  file.path(data_dir, "AdditionalFile26MethylationDataExample55.csv")
)
```

### 3. DNAm age 산출

RStudio에서 `scripts/r/Methyl_age_230526.R`를 열어 경로를 수정한 뒤 실행합니다. 완료되면 샘플별 DNAm age가 포함된 `Output.csv`와 검증 로그 `LogFile.txt`가 생성됩니다. 결과 파일은 프로젝트 구조를 유지하기 위해 `outputs/`에 저장하는 것을 권장합니다.

### 4. 노트북 분석

```bash
jupyter lab
```

다음 순서로 노트북을 실행합니다.

1. `notebooks/2_preprocessing.ipynb`
2. `notebooks/4_BDI_and_AgeDiff_Correlation_Analysis.ipynb`

각 노트북에 포함된 로컬 데이터 경로도 실행 환경에 맞게 변경해야 합니다.

## 결과 해석

`age_diff`는 일반적으로 DNAm age와 실제 연령의 차이로 사용됩니다. 계산 방향(`DNAm age - chronological age` 또는 그 반대)은 분석 전 명시하고 전체 분석에서 일관되게 유지하세요. 양수·음수의 해석은 이 정의에 따라 달라집니다.

BDI 분석은 상관계수, 정규성 검정, 등분산성 검정, 두 집단 비교를 포함합니다. 관찰된 관계는 표본 구성, 공변량, 데이터 품질의 영향을 받을 수 있으므로 인과관계로 해석하지 않아야 합니다.

## 데이터 및 재현성 관리

- 원천 데이터와 대용량 산출물은 Git에 커밋하지 않습니다.
- 개인식별 가능 정보나 민감한 임상 데이터는 저장소에 추가하지 않습니다.
- 입력 데이터 버전, 실행 날짜, 사용한 패키지 버전, 경로 변경 사항을 분석 기록에 남깁니다.
- 출력물은 `outputs/`에 저장하고, 재현 가능한 원본 코드와 분리해 관리합니다.

## 참고 자료

- Horvath, S. (2013). *DNA methylation age of human tissues and cell types*. Genome Biology, 14, R115.

## 개발 운영 방식

개인 프로젝트이지만 기능 브랜치, Issue, Pull Request를 기준으로 변경 이력을 관리합니다. 자세한 규칙은 [CONTRIBUTING.md](CONTRIBUTING.md)와 [GitHub 운영 가이드](docs/github-workflow.md)를 참고하세요.

## 외부 코호트 점검

후보 코호트의 조직·플랫폼·연령 메타데이터를 먼저 점검한 뒤에만 대용량 beta value를 내려받습니다. 후보와 제외 근거, 실행 방법은 [외부 코호트 선정과 메타데이터 점검](docs/cohort-selection.md)에 정리했습니다.

GSE40279의 beta value 수집 전 준비와 명시적 다운로드 방법은 [GSE40279 데이터 수집 준비](docs/gse40279-data-preparation.md)를 참고하세요.

과적합을 줄이기 위한 `GSE207605` 다중 코호트 선별 기준, 1차 학습 풀, 잠금 외부 검증 정책은 [다중 코호트 학습 풀과 과적합 방지 정책](docs/multicohort-training-pool.md)에 정리했습니다.

최종 Elastic Net 모델을 학습 데이터만으로 고정하고 `GSE87571`에서 외부 평가하는 순서와 재실행 방지 정책은 [최종 모델 고정과 외부 검증](docs/final-model-and-external-validation.md)을 참고하세요.

학습 표본 수가 성능의 병목인지 확인하는 코호트 단위 학습곡선 설계와 실행 결과는 [학습 데이터 규모별 학습곡선](docs/learning-curve-analysis.md)을 참고하세요.

고령층 과소예측과 취약 코호트의 연령 지원도·CpG 분포 이동 감사 결과는 [연령 극단과 취약 코호트 편향 감사](docs/age-cohort-bias-audit.md)를 참고하세요.

연령 밀도 가중치·혼합 가중치·내부 보정을 중첩 교차검증으로 비교한 결과와 새 외부 코호트 검증 전까지의 적용 범위는 [연령 편향 완화 후보 비교 실험](docs/age-bias-mitigation-experiment.md)을 참고하세요.
