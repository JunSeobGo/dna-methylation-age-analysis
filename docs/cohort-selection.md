# 외부 코호트 선정과 메타데이터 점검

## 목적

기존 `GSE72680`의 표본 부족 한계를 보완하기 위해, 혈액 조직과 Illumina HumanMethylation450K 플랫폼이 일치하는 공개 코호트를 우선 점검한다. 이 단계에서는 원천 beta value를 내려받지 않고 표본 메타데이터만 수집한다.

## 후보 코호트

| 코호트 | 제안 역할 | 선정 근거 |
| --- | --- | --- |
| [GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) | 학습 후보 | 전체혈액 450K, 넓은 연령 범위 |
| [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571) | 외부 검증 후보 | 전체혈액 450K, 14–94세 연령 범위 |
| [GSE50660](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50660) | 민감도 분석 후보 | 말초혈액 450K, 흡연 상태 점검 가능 |
| [GSE55763](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55763) | 규모 확장 후보 | 대규모 말초혈액 450K, 기술 반복 표본 확인 필요 |

## 기본 학습에서 제외하는 코호트

`GSE80417`, `GSE84727`, `GSE111629`은 혈액 450K 조건에는 맞지만 질환군과 대조군으로 구성되어 있다. 일반적인 연령 예측 모델의 기본 학습 자료에는 포함하지 않고, 향후 질환군 민감도 분석이 필요할 때 별도로 검토한다.

## 실행 방법

먼저 registry 형식만 확인하려면 다음을 실행한다.

```bash
Rscript scripts/r/2_cohort_metadata_audit.R --dry-run
```

실제 메타데이터를 수집하려면 `GEOquery`가 필요하다.

```r
install.packages("BiocManager")
BiocManager::install("GEOquery")
```

그다음 프로젝트 루트에서 실행한다.

```bash
Rscript scripts/r/2_cohort_metadata_audit.R
```

결과는 Git에서 제외되는 `outputs/cohort_metadata_audit/`에 저장된다.

## 통과 기준

- 실제 표본의 95% 이상에서 연령 정보를 확인할 수 있어야 한다.
- 실제 표본의 95% 이상이 registry의 450K 플랫폼·혈액 조직 조건과 일치해야 한다.
- 질환군 여부, 기술 반복 표본, 결측값 비율을 보고 수동 검토한다.
- 이 기준을 통과한 코호트만 beta value 다운로드와 정규화 단계로 이동한다.
