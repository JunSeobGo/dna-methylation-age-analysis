# DNA 메틸화 기반 연령 예측 프로젝트 개요

## 목적

여러 공개 혈액 DNA 메틸화 코호트의 beta value로 실제 연령을 예측하고, 새로운 코호트와 극단 연령대에서도 성능이 유지되는지 검증합니다. 단순한 평균 정확도뿐 아니라 코호트 이동, 연령 편향, 반복측정·가족 의존성과 데이터 누수를 함께 관리하는 것을 목표로 합니다.

본 프로젝트는 연구·교육용 분석이며 임상 진단이나 개인 건강 판단에 사용할 수 없습니다.

## 현재 상태

| 구분 | 상태 |
| --- | --- |
| 학습 데이터 | 8개 코호트, 5,541명, 공통 869 CpG |
| 내부 평가 | 코호트 그룹 기반 중첩 leave-one-cohort-out |
| 외부 검증 완료 모델 | 표준 Elastic Net, GSE87571 729명에서 평가 완료 |
| 개선 후보 | 연령 밀도 가중치 `age_density`, 학습 데이터로 잠금 완료 |
| 후보 외부 검증 | SATSA beta 추출 1회 실행, 사전 품질 게이트 `FAIL`(869 CpG 중 610개만 존재)로 평가 미수행 |
| 869 CpG 트랙 | SATSA 주 검증 수행 불가로 종결, 재시도하지 않음 |
| 610 CpG 제한 트랙 | SATSA 가용 특징만으로 모델 쌍을 재잠금해 재검증 예정 |
| 적용 경계 | 표준 Elastic Net을 현행 모델로 유지, `age_density`는 미검증 후보 |

세부 성능 수치는 [최종 모델과 외부 검증](final-model-and-external-validation.md) 및 [연령 편향 완화 실험](age-bias-mitigation-experiment.md)을 참고합니다.

## 전체 분석 흐름

```text
공개 코호트 조사
  → 메타데이터·질환·플랫폼 감사
  → include / conditional / exclude 정책 고정
  → 8개 코호트 공통 869 CpG 학습 bundle 생성
  → 코호트 그룹 기반 중첩 교차검증
  → 표준 모델 잠금
  → GSE87571 외부 검증 1회
  → 학습곡선과 연령·코호트 편향 감사
  → age-density 후보 비교·잠금
  → 새 SATSA 외부 검증 준비
  → SATSA beta 스트리밍 추출과 품질 게이트 판정(CpG coverage 미달로 FAIL)
  → 869 CpG 트랙 종결, 610 CpG 제한 트랙으로 재검증 개시
  → 사전 조건 통과 시에만 모델 교체 검토
```

SATSA 추출 결과와 실패 원인은 [SATSA beta 스트리밍 추출과 품질 게이트](satsa-beta-extraction.md)에 정리했습니다.

## 핵심 데이터 역할

| 데이터 | 역할 | 모델 선택 사용 |
| --- | --- | --- |
| GSE207605에서 선별한 8개 코호트 | 학습·내부 검증 | 사용 |
| GSE87571 | 표준 모델의 잠금 외부 검증 | 사용하지 않음 |
| SATSA S-BSST1206 | age-density 후보의 새 고령층 외부 검증 | 사용하지 않음 |
| GSE193879 건강 대조군 | 향후 소아 민감도 분석 후보 | 사용하지 않음 |
| conditional 코호트 | 반복·질환 조건을 통제한 민감도 분석 | 기본 모델에는 사용하지 않음 |

코호트별 포함·제외 근거는 `config/`의 CSV manifest와 [다중 코호트 학습 풀](multicohort-training-pool.md)에서 관리합니다.

## 재현성과 누수 방지 원칙

1. 행 단위 무작위 분할 대신 원 코호트 단위로 학습과 평가를 분리합니다.
2. 결측 대치, 표준화와 하이퍼파라미터 선택은 각 학습 fold 안에서만 수행합니다.
3. 외부 데이터는 모델과 판정 기준을 잠근 뒤에만 읽습니다.
4. 한 번 평가에 사용한 외부 코호트를 새 후보 선택에 재사용하지 않습니다.
5. 전체 MAE와 함께 연령 구간·코호트별 오차, 편향과 calibration을 보고합니다.
6. 반복측정과 가족 표본은 독립 표본으로 취급하지 않습니다.
7. 원본·처리 데이터와 생성 결과는 Git에서 제외하고 URL·크기·해시와 코드를 기록합니다.

결정 배경과 대안은 [ADR 인덱스](adr/README.md)에서 확인할 수 있습니다.

## 저장소 구조

```text
DNAage_clean/
├── config/                 # 코호트 정책과 다운로드 manifest
├── data/
│   ├── raw/                # 공개 원본, Git 제외
│   └── processed/          # 분석 bundle·잠금 모델, Git 제외
├── docs/
│   ├── adr/                # 기술·분석 의사결정 기록
│   ├── figures/            # 저장소에 포함하는 결과 그림
│   └── README.md           # 문서 인덱스
├── notebooks/              # 초기 Python 탐색 분석
├── outputs/                # 재생성 가능한 지표·그림, Git 제외
├── scripts/
│   ├── ci/                 # 저장소 정책 검사
│   └── r/                  # 번호 순서의 현행 R 파이프라인과 legacy/
└── .github/                # Issue·PR 양식과 CI
```

## 문서 진입점

- [전체 문서 지도](README.md)
- [R 스크립트 실행 순서](../scripts/r/README.md)
- [GitHub 작업 방식](github-workflow.md)
- [기여 가이드](../CONTRIBUTING.md)
- [ADR 인덱스](adr/README.md)
