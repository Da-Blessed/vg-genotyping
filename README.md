# vg short-read and long-read genotyping workflow

Phased panel VCF로 variation graph를 만들고, short reads 또는 long reads를 Giraffe로 mapping한 뒤 panel allele의 genotype을 계산하는 Snakemake workflow입니다.

실제 실행 명령은 `vg call -v panel.vcf.gz`이지만 이 workflow는 **novel variant calling이 아니라 VCF-guided genotyping**입니다. `-v`로 panel VCF를 제공하여 이미 정의된 variant site와 allele을 대상으로 genotype을 산출합니다.

## 지원 read

- short-read: paired-end FASTA/FASTQ, Giraffe `default` preset
- PacBio HiFi long-read: Giraffe `hifi` preset
- Oxford Nanopore R10 long-read: Giraffe `r10` preset

long-read Giraffe 지원은 비교적 최근에 추가되어 활발히 개발 중입니다. vg 버전과 preset을 고정하고 실제 dataset으로 정확도와 resource를 검증하는 것을 권장합니다. 이 workflow는 검증한 vg `1.70.0`의 long-read index 파일명과 CLI를 사용합니다.

## Panel genotype 검증

기본값인 `validation.strict: true`에서는 panel의 호출된 diploid genotype이
phased인지, allele index가 유효한 숫자인지를 검사합니다. Missing genotype이
있는 panel을 사용하되 나머지 strict 검사를 유지하려면 다음과 같이 설정합니다.

```yaml
validation:
  strict: true
  allow_missing_genotypes: true
```

이 설정은 `.|.`, `./.` 및 `0|.` 같은 missing genotype을 허용하지만, `0/1`
같은 unphased called genotype은 계속 거부합니다. `allow_missing_genotypes`의
기본값은 `false`입니다. VCF 정렬 순서, reference 범위, record 간 중첩 및
sequence-resolved allele 검사는 이 옵션과 관계없이 항상 수행됩니다.

## 입력

- 압축하지 않은 reference FASTA
- 압축하지 않은 phased, sequence-resolved multi-sample panel VCF
- short-read sample sheet
- long-read sample sheet

```yaml
reference: /path/to/reference.fa
panel_vcf: /path/to/panel.vcf
shortread_samples: config/samples.shortread.tsv
longread_samples: config/samples.longread.tsv
```

Short-read sheet:

```tsv
sample	read1	read2
sampleA	/path/to/A_R1.fastq.gz	/path/to/A_R2.fastq.gz
```

Long-read sheet의 `platform`은 `hifi` 또는 `r10`입니다.

```tsv
sample	read	platform
sampleB	/path/to/B.hifi.fastq.gz	hifi
sampleC	/path/to/C.ont.fastq.gz	r10
```

## 설치와 실행

```bash
conda env create -f environment.yaml
conda activate vg-genotyping-workflow

# short-read와 long-read 모두
snakemake --profile profiles/default

# short-read만
snakemake --snakefile workflow/shortread.smk --profile profiles/default

# long-read만
snakemake --snakefile workflow/longread.smk --profile profiles/default
```

먼저 `--dry-run`을 붙여 DAG와 입력 경로를 확인할 수 있습니다.

## 처리 구조

1. reference와 panel VCF 검증
2. `vg construct -a`로 VCF allele을 보존한 graph 생성
3. XG, GBWT, GBZ 및 read-type별 Giraffe index 생성
4. Giraffe mapping과 `vg pack`
5. XG graph에 `vg call -v panel.vcf.gz`를 적용해 sample별 genotyping
6. sample VCF를 cohort VCF로 병합

## 결과

| 경로 | 내용 |
|---|---|
| `results/shortread/vcf/{sample}.vcf.gz` | short-read sample genotype VCF |
| `results/shortread/cohort/vg.shortread.cohort.vcf.gz` | short-read cohort VCF |
| `results/longread/vcf/{sample}.vcf.gz` | long-read sample genotype VCF |
| `results/longread/cohort/vg.longread.cohort.vcf.gz` | long-read cohort VCF |

SnpEff와 AnnotSV annotation은 이 저장소에 결합하지 않았습니다. 생성한 cohort VCF를 독립 `variant-annotation` workflow의 `config/inputs.tsv`에 등록해 실행합니다.

## 테스트

```bash
python -m unittest discover -s tests -v
snakemake --snakefile Snakefile --configfile tests/config.yaml --lint
snakemake --snakefile Snakefile --configfile tests/config.yaml --cores 1 --dry-run
snakemake --snakefile workflow/shortread.smk --configfile tests/config.yaml --cores 1 --dry-run
snakemake --snakefile workflow/longread.smk --configfile tests/config.yaml --cores 1 --dry-run
```

참고 문서: [vg automatic indexing](https://github.com/vgteam/vg/wiki/Automatic-indexing-for-read-mapping-and-downstream-inference), [long-read Giraffe](https://github.com/vgteam/vg/wiki/Mapping-long-reads-with-Giraffe), [SV genotyping and variant calling](https://github.com/vgteam/vg/wiki/SV-Genotyping-and-variant-calling)
