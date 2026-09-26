# Beginner and interview guide

## The shortest correct explanation

This project helps a reader move from a **type 2 diabetes GWAS signal** to a more organised view of the variants, genes, and molecular evidence around that signal.

A GWAS usually identifies a region of DNA associated with a disease. It does not automatically identify the responsible variant or gene. This explorer brings together three existing Open Targets evidence layers for five selected FinnGen regions:

1. variants retained by fine-mapping;
2. genes ranked by the Open Targets Locus-to-Gene model; and
3. molecular-QTL colocalisation records.

The project checks and organises these public records. It does not diagnose diabetes, predict whether a person will develop diabetes, or prove that a gene causes diabetes.

## Is this regression or classification?

### Direct answer

**The explorer itself is neither a regression model nor a classification model.** It is a deterministic **post-GWAS evidence-integration and audit pipeline** followed by a descriptive case study and interactive visualisation.

There is no target variable that this repository learns to predict. It does not train a model, split data into training and test sets, or calculate predictive accuracy. Its own work is to retrieve public aggregate results, verify the saved responses, select five regions using a documented rule, calculate genomic distances, match records by gene identifiers, and display existing source scores without combining them.

### Why regression and classification still appear in the background

Different upstream stages use different methods:

| Stage | Question | Method type | Did this repository run it? |
|---|---|---|---|
| FinnGen type 2 diabetes GWAS | Is each genetic variant associated with type 2 diabetes case/control status? | Binary-trait association testing using logistic mixed-model methods | No |
| Open Targets fine-mapping | Which variants could account for an association signal? | Probabilistic statistical inference | No |
| Open Targets L2G | Which nearby gene should be prioritised for a credible set? | Supervised gradient-boosting classification during training; the resulting score is used for gene ranking | No |
| Open Targets colocalisation | Could a disease signal and a molecular-QTL signal share an underlying variant? | Probabilistic hypothesis comparison | No |
| This explorer | What evidence was returned, how do the records match, and what is missing? | Data retrieval, validation, evidence integration, descriptive analysis, and visualisation | Yes |

The important interview distinction is that **logistic regression is the statistical machinery behind a binary-trait GWAS, but the goal is association testing rather than building a patient classifier**. The output asks whether variant dosage is associated with disease status after accounting for relevant structure and covariates. It is not a tool that receives one new patient's information and predicts that patient's diagnosis.

## The problem statement

### In everyday language

A GWAS can tell us that a neighbourhood of the genome is related to type 2 diabetes. That neighbourhood may contain many correlated variants and several genes. The difficult question is: **which variant and gene should researchers investigate next, and how strong is the supporting evidence?**

No single displayed number answers that question. Fine-mapping, L2G, genomic distance, and molecular-QTL colocalisation describe different parts of the evidence. Missing records also cannot safely be changed into zeros. The project therefore organises the evidence, preserves its source and limitations, and lets a reader compare the layers without presenting them as proof of causality.

### In data-science language

**Input:** Public aggregate Open Targets records for five selected credible sets from `FINNGEN_R12_T2D`, including variant PIPs, returned L2G candidates, and molecular-QTL colocalisation rows.

**Processing:** Response-hash verification, pagination and count checks, deterministic locus selection, probability validation, lead-variant-to-TSS distance calculation, and candidate-to-QTL matching with Ensembl gene identifiers.

**Output:** Validated source-derived tables, an audit report, a five-region descriptive case study, and an interactive explorer.

**Research question:** For five selected FinnGen type 2 diabetes regions, what variant, candidate-gene, and molecular-QTL evidence did Open Targets return, where do those records agree, and where are the evidence gaps?

**Non-goals:** Patient-level prediction, causal-gene confirmation, a new GWAS, a new L2G model, clinical target ranking, or experimental validation.

## Start with the biological foundation

### DNA

DNA is the molecule that stores biological instructions. Its sequence is written with four chemical bases, represented as **A, C, G, and T**. The order of these bases carries information.

**Simple picture:** DNA is the text used to write the body's instruction manual.

### Genome

A genome is the complete set of DNA instructions in a person or organism. In humans, most of this DNA is organised into 23 pairs of chromosomes in the cell nucleus, with a small additional genome in mitochondria.

**Simple picture:** If DNA is text, the genome is the complete collection of volumes.

### Chromosome

A chromosome is one long DNA molecule packaged with proteins. Chromosome numbers help identify where a variant is located. A variant written as `10_114758349_C_T`, for example, is on chromosome 10 at a particular genomic position.

**Simple picture:** A chromosome is one volume in the complete genome collection.

### Gene

A gene is a region of DNA whose information is used to make a functional RNA or protein. Proteins and RNAs help cells perform biological work. A disease-associated variant can be inside a gene, near a gene, or in a regulatory region that affects a more distant gene.

**Do not say:** “The closest gene must be the causal gene.” Physical distance is one clue, not proof.

### Gene expression

Gene expression describes how the information in a gene is used to produce RNA. Researchers often measure RNA abundance as an indication of how active a gene is in a particular tissue or cell type.

**Simple picture:** The gene is an instruction; expression describes how much that instruction is being read.

### Variant and SNP

A genomic variant is a DNA sequence difference among people. A **single-nucleotide polymorphism**, or **SNP** (pronounced “snip”), is a difference at one DNA base position.

A variant can be associated with disease without directly causing it. It may simply be inherited together with the causal variant.

### Linkage disequilibrium

Nearby variants are often inherited together. This correlation is called **linkage disequilibrium**, usually shortened to **LD**. Because of LD, several variants in one region can show a similar disease association even if only one, or possibly none of the measured variants, directly changes the biological mechanism.

**Why it matters here:** LD is one reason a GWAS signal initially points to a region rather than one certain variant.

### Locus

A locus is a position or neighbourhood in the genome. In this project, a locus means the region surrounding a GWAS association signal and its credible set.

**Simple picture:** The genome gives the complete map; a locus is one marked neighbourhood on that map.

### Phenotype

A phenotype is a measured trait or condition. It can be binary, such as type 2 diabetes case versus control, or quantitative, such as height, cholesterol, or gene-expression level.

## From a GWAS signal to plausible variants

### GWAS

GWAS means **genome-wide association study**. A GWAS examines variants across the genomes of many people and tests whether each variant is statistically associated with a trait.

For a binary disease such as type 2 diabetes, the analysis compares case and control status while accounting for covariates and population or family structure. A small p-value indicates that the observed association would be unusual under the model's no-association assumption. It does not prove that the variant causes the disease.

### Association is not causation

An association means that a variant and a trait vary together in the analysed population. The variant may be causal, may tag another nearby causal variant through LD, or may be affected by study-specific limitations. Establishing a biological mechanism requires more evidence.

### P-value

A p-value measures how incompatible the observed data are with a specified null model, usually a model of no association. A very small p-value supports evidence of an association signal.

**It is not:**

- the probability that the variant is causal;
- the probability that the finding is “just chance”; or
- the size or clinical importance of the effect.

### Lead variant

The lead variant is the representative variant for an association signal, often the one with the strongest reported evidence or highest fine-mapping support under the source definition. It is a useful label for the region. It is not automatically the causal variant.

### Fine-mapping

Fine-mapping is a statistical analysis performed after a GWAS signal is identified. It uses association evidence and the correlation among variants to narrow a broad signal to a smaller group of plausible variants.

**Simple picture:** GWAS identifies the neighbourhood; fine-mapping produces a shortlist of addresses within it.

This repository displays Open Targets fine-mapping outputs. It does not run its own fine-mapping model.

### PIP

PIP means **posterior inclusion probability**. In this context, it is the source fine-mapping model's support for including a variant as causal for that signal, given the data, model, and assumptions.

A higher PIP means the source model places more support on that variant relative to other variants in the credible set. PIP is not a probability that a gene causes type 2 diabetes, and it is not universal certainty outside the model.

### Credible set

A credible set is a group of fine-mapped variants whose combined PIP reaches a chosen level, commonly 95%. Under the fine-mapping model and assumptions, the set is intended to contain the causal variant with that level of probability.

A credible set with one high-PIP variant is concentrated. A set with many low-PIP variants is less resolved because several variants remain plausible.

## From plausible variants to candidate genes

### TSS and genomic distance

TSS means **transcription start site**, the genomic position where transcription of a gene begins. The explorer calculates the absolute distance between a lead variant and the canonical TSS of each returned candidate gene.

Distance can help prioritise genes, but a nearby variant can regulate a distant gene. The nearest returned gene is therefore a comparison, not a causal conclusion.

### Ensembl gene identifier

An Ensembl gene identifier is a stable database identifier such as `ENSG00000148737`. Gene symbols can change or be ambiguous, so the explorer uses these identifiers when matching candidate genes to molecular-QTL records.

### Locus-to-Gene and the L2G score

Open Targets' **Locus-to-Gene**, or **L2G**, model ranks protein-coding genes near a GWAS credible set. It uses features such as genomic distance, molecular-QTL colocalisation, and predicted variant consequences.

The model is trained with positive and negative gene–credible-set examples using gradient boosting. Its 0-to-1 L2G output is used to rank candidate genes. The explorer retrieves and displays that existing score; it did not train or validate the L2G model.

**Do not say:** “An L2G score of 0.8 proves an 80% chance that this gene causes diabetes.” It is safer to say that Open Targets gave the gene a higher model score for that credible set.

## QTL and molecular evidence

### Quantitative trait

A quantitative trait is something measured on a numerical scale, such as gene-expression level, protein abundance, or the proportion of a gene's RNA that is spliced in a particular way.

### QTL

QTL means **quantitative trait locus**. It is a genomic location where genetic variation is statistically associated with variation in a measured quantitative trait.

Suppose people with one version of a variant tend to produce more RNA from a gene than people with another version. That variant may be part of an expression QTL for the gene.

QTL does not mean “a diabetes gene.” It means that a genetic signal is associated with a measured molecular trait.

### Common molecular-QTL types

- **eQTL:** associated with gene-expression level.
- **pQTL:** associated with protein abundance.
- **sQTL:** associated with RNA splicing.
- **tuQTL:** associated with transcript usage.

These effects can depend on tissue or cell type. An eQTL found in blood may not behave the same way in liver, pancreas, muscle, or another relevant tissue.

### Colocalisation

Colocalisation compares two association signals, such as a type 2 diabetes GWAS signal and an eQTL signal. It asks whether the patterns are consistent with a shared underlying variant.

**Simple picture:** Two alarms are sounding in the same genomic neighbourhood. Colocalisation asks whether one source may be triggering both alarms or whether separate sources are more likely.

Colocalisation is stronger evidence than simply noticing that two signals are nearby, but it does not prove that changing the gene causes or prevents disease.

### H3 and H4

For the COLOC-PIP output used by Open Targets:

- **H3** supports a model in which both traits have association signals but different variants account for them.
- **H4** supports a model in which the two traits share an underlying variant.

Higher H4 is more consistent with a shared-variant model in that analysis. It is not a universal causal probability and must be interpreted with tissue, study, fine-mapping, and model assumptions.

### CLPP

CLPP is the colocalisation posterior probability statistic produced by the eCAVIAR approach. It is based on variant-level fine-mapping probabilities from the two signals.

CLPP and H4 come from different methods. The explorer keeps them in separate columns and does not add or average them.

### Candidate-matched and unmatched QTL rows

The explorer calls a colocalisation row **candidate-matched** only when its source gene identifier matches an Ensembl identifier in the returned L2G candidate list. Other rows remain visible as locus-level evidence.

This prevents a large number of QTL rows at a locus from being incorrectly described as support for one particular candidate gene.

## Data and engineering terms

### Aggregate data

Aggregate data summarise study results without providing individual participant records. This project uses public aggregate Open Targets records. It does not contain personal medical histories or participant-level genotypes.

### API

An application programming interface, or API, is a structured way for software to request data from another system. The project used the Open Targets GraphQL API to retrieve the saved source records.

### Pagination

An API may return a large result in several pages. The project checks that all pages advertised for each bounded query were retrieved and that rows were not duplicated.

### Hash

A SHA-256 hash is a digital fingerprint of a file. If a cached response changes by even a small amount, its hash changes. The pipeline verifies the saved hashes before rebuilding the explorer.

### Provenance

Provenance records where data came from and how it was processed. The project saves each request, variables, retrieval time, source URL, byte count, and response hash so displayed evidence can be traced back to the saved response.

### Missing, not returned, and zero

These are different:

- **Zero:** the source measured or calculated a value of zero.
- **Missing:** a field has no value.
- **Not returned:** the bounded source response did not include a corresponding record.

The explorer preserves these differences. It does not turn missing evidence into negative evidence.

## What the project actually did

1. Retrieved the study metadata and 368 available credible sets for the FinnGen R12 type 2 diabetes study from Open Targets release 26.09.
2. Selected five regions deterministically using ascending source p-value, stable tie-breaking, and at least one megabase between same-chromosome lead variants.
3. Retrieved the source variants, L2G candidates, and molecular-QTL colocalisation records for those regions.
4. Verified nine cached API responses, response hashes, pagination, identifiers, counts, probability ranges, and PIP sums.
5. Calculated lead-variant-to-TSS distances for returned candidates when coordinates were compatible.
6. Matched QTL records to returned candidates using Ensembl gene identifiers while retaining unmatched locus-level records.
7. Produced 65 variant rows, 11 candidate-gene rows, and 44 molecular-QTL colocalisation rows.
8. Kept PIP, L2G, H4, and CLPP separate and preserved missing values and source quality-control warnings.
9. Generated downloadable tables, an audit report, a case study, and the interactive explorer.

## What to say in an interview

### Thirty-second version

> A GWAS can identify a genomic region associated with type 2 diabetes, but it usually cannot directly tell us which variant or gene is responsible. I built a post-GWAS evidence explorer that brings together Open Targets fine-mapped variants, Locus-to-Gene rankings, and molecular-QTL colocalisation for five FinnGen regions. My pipeline verifies the saved source data, matches molecular records to candidate genes, preserves missing values and warnings, and keeps the different statistical scores separate. It is a descriptive evidence-integration project, not a patient classifier or a causal-gene model.

### If asked whether it is regression or classification

> My repository itself is neither. It does not train a predictive model. It retrieves, validates, matches, and presents existing post-GWAS evidence. The upstream type 2 diabetes GWAS is a binary-trait association analysis that uses logistic mixed-model methods. Open Targets' upstream L2G model is trained as a gradient-boosting classification model and its score is used to rank candidate genes. Fine-mapping and colocalisation are probabilistic inference steps. I kept those upstream methods distinct from what my own pipeline implemented.

### If asked what a QTL is

> A QTL is a genomic location where a genetic variant is associated with a numerical biological measurement. For example, an eQTL links genetic variation to the amount of RNA produced by a gene. In this project, colocalisation checks whether a diabetes GWAS signal and a molecular-QTL signal may share an underlying variant. That can support a biological hypothesis, but it does not prove that the gene causes diabetes.

## Claims to avoid

- “I trained a diabetes classifier.”
- “I ran the FinnGen GWAS.”
- “I created the L2G model.”
- “The top gene is the causal gene.”
- “H4, CLPP, PIP, and L2G are all the same kind of probability.”
- “No returned QTL record means there is no molecular relationship.”
- “The five selected regions represent all cardiometabolic diseases.”
- “These data contain individual FinnGen patients.”

## Primary references

- [NHGRI genome and genetics glossary](https://www.genome.gov/genetics-glossary)
- [NHGRI explanation of GWAS](https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies-GWAS)
- [FinnGen guidance on binary and quantitative GWAS models](https://docs.finngen.fi/working-in-the-sandbox/running-analyses-in-sandbox/how-to-run-genome-wide-association-studies-gwas)
- [Open Targets credible sets and variant PIP](https://platform-docs.opentargets.org/credible-set)
- [Open Targets Locus-to-Gene methodology](https://platform-docs.opentargets.org/gentropy/locus-to-gene-l2g)
- [Open Targets colocalisation methods](https://platform-docs.opentargets.org/gentropy/colocalisation)
