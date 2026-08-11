# GitHub 운영 가이드

## 운영 원칙

하나의 Issue와 Pull Request에는 하나의 의사결정과 검증 가능한 결과만 남깁니다. 개인 프로젝트라도 작업 근거와 변경 이력을 추적할 수 있어 포트폴리오 설명에 유리합니다.

## Issue

각 Issue에는 목적, 작업 목록, 검증 가능한 완료 조건을 작성합니다. 구현 Pull Request에서는 `Closes #<Issue 번호>`로 연결합니다.

| 라벨 | 용도 |
| --- | --- |
| `type: feature` | 새 파이프라인 기능 또는 분석 단계 |
| `type: bug` | 재현성, 계산, 문서 오류 |
| `type: documentation` | README 또는 분석 문서 개선 |
| `type: refactor` | 의도한 동작을 바꾸지 않는 구조 개선 |
| `status: ready` | 작업을 시작할 수 있는 상태 |
| `priority: high` | 다음 작업 주기에 우선 처리할 항목 |

## Pull Request

- 기능·수정 Pull Request는 `develop`을 대상으로 생성합니다.
- 릴리스 Pull Request는 `develop`에서 `main`으로 생성합니다.
- Pull Request에 목적, 변경 사항, 검증 결과, 관련 Issue를 작성합니다.
- 병합 뒤 작업 브랜치는 삭제합니다.

Pull Request에서는 저장소 파일·설정 CSV 검증과 R 스크립트 구문 검사가 각각 독립 작업으로 실행됩니다. 이 검사는 원천 데이터나 분석 패키지를 요구하지 않으므로 데이터 다운로드 여부와 관계없이 빠르게 실패 원인을 확인할 수 있습니다.

## Milestone과 Release

| Milestone | 완료 기준 |
| --- | --- |
| `v0.1` | 재현 가능한 DNAm age 기준 파이프라인 |
| `v0.2` | 다중 코호트 품질검사와 외부 검증 |
| `v1.0` | 분석 보고서와 재현 가능한 배포 버전 |

각 milestone이 끝나면 `main`에 병합하고 버전 태그와 GitHub Release를 생성합니다.

## GitHub Projects

활성 Issue가 다섯 개 이상 생기면 Project 보드를 만듭니다. `Todo → In Progress → Review → Done` 흐름을 사용하고, 별도 할 일 목록 대신 Issue를 연결해 관리합니다.
