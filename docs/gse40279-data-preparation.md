# GSE40279 데이터 수집 준비

## 범위

GSE40279는 656개 전체혈액 표본의 Illumina HumanMethylation450K 데이터다. 공개된 평균 beta value 파일은 압축 상태로 약 1.1GB이므로, 메타데이터 감사와 달리 명시적인 다운로드 명령에서만 수집한다.

## 파일 구성

| 리소스 | 크기 | 용도 |
| --- | ---: | --- |
| `GSE40279_sample_key.txt.gz` | 약 4.5KB | GSM 접근번호와 beta matrix 열 이름 연결 |
| `GSE40279_average_beta.txt.gz` | 약 1.1GB | 656개 표본의 평균 beta value 행렬 |

원격 URL과 로컬 저장 위치는 [`config/dataset_download_manifest.csv`](../config/dataset_download_manifest.csv)에서 관리한다. `data/raw/`는 Git에서 제외된다.

## 실행 방법

기본 실행은 파일을 내려받지 않고 다운로드 계획만 표시한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/3_prepare_gse40279_data.R
```

먼저 작은 표본 키 파일만 내려받는다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/3_prepare_gse40279_data.R --download-sample-key
```

저장 공간을 확인한 뒤 평균 beta value 행렬을 받는다. 이 명령은 표본 키도 함께 확인한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/3_prepare_gse40279_data.R --download-beta
```

다운로드는 `.part` 임시 파일에 수행하고 gzip 형식을 확인한 뒤 최종 파일명으로 바꾼다. 네트워크가 끊기면 `.part` 파일을 유지하며, 같은 명령을 다시 실행하면 중단된 지점부터 재개한다. 처음부터 다시 받으려면 `--force`를 추가한다.

## 표본 매핑

표본 키의 가운데 열 participant ID는 GEO 표본 title의 마지막 숫자와 대응한다. beta matrix 열은 participant ID 앞에 `X`가 붙은 형식(예: `1001` → `X1001`)이다. 표본 키의 마지막 열은 별도 배열 칩 ID로 보존한다. 다음 명령으로 beta column 이름·GSM 접근번호·실제 연령을 하나의 매니페스트로 만든다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/4_build_gse40279_sample_manifest.R
```

결과 `data/processed/GSE40279_sample_manifest.csv`에는 656개 표본의 `beta_column_name`, GSM 접근번호, 실제 연령이 포함된다. 이 파일은 beta matrix 열 검증과 모델 학습 라벨로 사용한다.

## 다음 검증

1. beta matrix의 첫 열이 CpG probe ID인지 확인한다.
2. 나머지 열이 `beta_column_name`과 정확히 일치하는지 점검한다.
3. 656개 표본 전체가 실제 연령과 연결된 뒤에만 학습용 행렬을 만든다.
