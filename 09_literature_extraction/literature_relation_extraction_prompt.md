# PlantQTLdb Literature-Based Biological Relation Extraction Prompt

**Model:** DeepSeek V4 Pro

**Extraction conventions:** Species-aware gene symbols and single-gene subjects.

> Reader note: The prompt below refers to paper text, paper identifiers, species metadata and a JSON output schema supplied by the caller. The caller-supplied schema is not included in this file. This note is documentation, not part of the extraction instructions.

## Core principle

Work only from the provided text. Do not use outside knowledge.
If the text does not clearly support a relation, do not extract it.

Your goal is high-recall but evidence-controlled extraction for PlantQTLdb.
Extract biologically meaningful gene-centered relations from the paper, including
expression, localization, trait/phenotype, stress response, downstream regulation,
physiological readouts, no-effect mutant results, and gene-family information.


---

## Compact metadata rule for large gene-family papers

For large gene-family papers, save output budget for the actual `relations` list.

- `GENE_BELONGS_TO_FAMILY` relations should still be exhaustive when the paper clearly
  defines a full member range or full member list.
- Do NOT exhaustively enumerate every family member in `detected_entities`.
  Use a compact range or summary there, e.g. `ONAC001-ONAC075`.
- Do NOT exhaustively enumerate every family member in `gene_id_links`.
  Include gene_id_links only for experimentally characterized focal genes or a few
  explicitly discussed representative genes. For large genome-wide family tables,
  leave gene_id_links empty or concise.
- Do not duplicate the same long family list across detected_entities, gene_id_links,
  and relations. The exhaustive family membership belongs in relations only.

---

## Output discipline

Use ONLY the JSON schema requested by the caller.
Return valid JSON only, no markdown and no commentary.
Every requested paper_id must appear exactly once.
If a paper has no extractable relation, return an empty relations list and a filled document_assessment.

Each relation must include a direct evidence_text sentence or sentence fragment from the provided document.
If a field is unknown, use null or an empty list.

---

## Gene symbol canonicalization discipline

Use species-aware canonical gene symbols, but do not invent names.

- Prefer canonical species-specific gene symbols when the species context is explicit
  and the paper itself supports the equivalence.
- For rice / Oryza sativa, prefer Os-prefixed symbols in `subject_gene` when the
  paper uses both an Os-prefixed form and a bare form for the same gene.
  Example: if the same rice paper uses both `OsDREB1` and `DREB1` for the same
  gene, use `OsDREB1` consistently in `subject_gene`; do not emit separate
  `OsDREB1` and `DREB1` relations for the same evidence.
- Do not add a species prefix just because the caller species is known. The paper
  must explicitly use or support the prefixed symbol.
- Do not merge orthologs across species. In comparative or multi-species papers,
  keep `OsDREB1A`, `AtDREB1A`, `ZmDREB1`, etc. separate unless the relation is
  explicitly about an ortholog group.
- Do not treat family names or subfamily names as individual genes.

Single-gene subject rule:
- `subject_gene` must normally contain exactly one gene symbol.
- Do not use compound subjects such as `OsMADS7 and OsMADS8` or
  `OsMADS1, OsMADS5, OsMADS7 and OsMADS8`.
- If the paper reports a double mutant, simultaneous knockdown, combined
  perturbation, or multi-gene construct, split the result into one relation per
  named gene and preserve the combination context in `condition_context`, e.g.
  `double mutant involving OsMADS7 and OsMADS8` or `simultaneous knockdown of
  OsMADS1, OsMADS5, OsMADS7 and OsMADS8`.
- If a statement is only about a gene family or broad unnamed group and no
  individual genes are named, use the family/group as the target or context, not
  as a compound `subject_gene`.

---

## Allowed relation_family values

Use ONLY these relation_family values. Do not invent labels such as "induced_by",
"upregulates", "localizes_to", "overexpression_phenotype", or "tissue_expression".
Put the original wording in raw_relation_wording instead.
Do NOT emit GENE_FAMILY_HAS_SIZE or GENE_FAMILY_IN_SPECIES. Family size can be derived later by counting GENE_BELONGS_TO_FAMILY relations; species is stored in the caller-provided species field.

1. GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE
2. GENE_BIOLOGICAL_PROCESS_OR_FUNCTION
3. GENE_RESPONDS_TO_CONDITION
4. GENE_EXPRESSED_IN
5. GENE_ASSOCIATED_WITH_METABOLITE_OR_PATHWAY
6. GENE_INTERACTS_WITH_GENE_OR_PROTEIN
7. GENE_ASSOCIATED_WITH_VARIANT_QTL_LOCUS
8. GENE_BELONGS_TO_FAMILY

---

## Focal-gene coverage checklist

For each focal gene in a paper, actively check whether the text supports each of
these categories. Extract them when directly supported:

1. Stress/hormone/treatment-induced expression or response:
   drought, salt/high salinity, cold/low temperature, heat, ABA, auxin,
   gibberellin, cytokinin, pathogen, submergence, nitrogen deficiency, light/dark,
   chemical treatment, inhibitor treatment.

2. Static tissue/cell-type expression:
   root, leaf, shoot, stem, panicle, spikelet, anther, pollen, seed, embryo,
   endosperm, coleoptile, seminal root, vascular tissue, phloem, xylem,
   meristem, epidermis, guard cell, tapetum, etc.

3. Subcellular localization:
   nucleus, cytoplasm, chloroplast, mitochondrion, membrane, endoplasmic reticulum,
   Golgi, plastid, cell wall, apoplast, etc.

4. Transgenic/mutant phenotype:
   overexpression, RNAi, CRISPR/Cas9, T-DNA insertion, knockout, loss-of-function,
   gain-of-function, complementation, natural allele, near-isogenic line.

5. Physiological or quantitative readouts directly supporting a gene phenotype:
   survival rate, photochemical efficiency, chlorophyll fluorescence, JIP test,
   photosynthetic rate, stomatal conductance, water loss rate, relative water content,
   electrolyte leakage, ROS, MDA, proline, ABA content, enzyme activity,
   plant height, tiller number, grain yield, grain size, grain width, grain length,
   fertility, flowering time.

6. Downstream genes, proteins, pathways, and transcriptome effects:
   RNA-seq, microarray, qRT-PCR validation, ChIP, EMSA, yeast one-hybrid,
   dual-luciferase, co-expression, protein interaction, TF-target regulation.

7. Negative/no-effect findings:
   no significant phenotype, no obvious difference, knockout had no effect,
   expression was not changed, not induced by treatment.
   Extract these when they answer a biologically meaningful gene-trait,
   gene-condition, or gene-function question.

Do not force all categories to appear. Extract only what is supported.

---

## Target relation families

### 1. GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE

Named agronomic trait, phenotype, physiological readout, or observed outcome.
Include: grain weight, chalkiness, amylose content, gel consistency, dwarfism,
drought tolerance, salt sensitivity, lodging resistance, flowering time, plant height,
survival rate, photochemical efficiency, water loss, ROS accumulation, grain yield.

Use this for mutant/transgenic phenotypes and direct physiological assay readouts
when tied to a gene perturbation.

Examples:
- "OsERF83 overexpression increased drought tolerance" →
  OsERF83 --[GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE]--> drought tolerance,
  experiment_basis="overexpression; drought_phenotyping", polarity="positive".
- "OsERF83OX plants showed reduced grain yield" →
  OsERF83 --[GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE]--> grain yield,
  experiment_basis="overexpression; yield_component_analysis", polarity="negative".
- "os erf83 mutants showed no obvious drought phenotype" →
  OsERF83 --[GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE]--> drought tolerance,
  experiment_basis="crispr_cas9_editing; drought_phenotyping", polarity="neutral",
  conclusion_text="knockout showed no obvious effect on drought tolerance".

Phenotype splitting rule — extract BOTH as separate relations only when clearly distinct:
- numerical measurement/readout: e.g., "AC as low as 0.6%", "survival rate 83–95%";
- named phenotypic class: e.g., "glutinous phenotype", "drought tolerance";
- statistical difference vs wild type belongs in evidence_text, not as its own relation.

Do NOT merge a numerical readout with a named phenotype into one relation when both
are independently meaningful and separately stated.

### 2. GENE_BIOLOGICAL_PROCESS_OR_FUNCTION

Gene/protein function, molecular activity, regulatory role, or pathway membership.
Include: transcription regulation, signal transduction, metabolic pathway function,
enzyme activity, kinase activity, protein complex formation, chromatin remodeling,
photosynthesis-related function, stress-response regulation.

Use this when the paper states a gene function/process, not just a phenotype.
For physiological readouts, prefer GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE unless
the text explicitly states a biological function.

### 3. GENE_RESPONDS_TO_CONDITION

ONLY when a treatment, stress, or condition changes gene expression, abundance,
activity, localization, or behavior.
The diagnostic question: "Does the text describe how a CONDITION affects this gene?"

Include: drought, high salinity/salt, ABA, cold, heat, nitrogen deficiency,
pathogen, submergence, hormone treatment, dark/light transition, UV irradiation,
elicitor treatment, inhibitor treatment.

Polarity rule:
- upregulated / induced / increased expression = positive
- downregulated / repressed / reduced expression = negative
- no induction / no change = neutral
- unclear wording = unclear

Examples:
- "OsERF83 was induced by drought, salt, cold and ABA" → emit separate relations:
  OsERF83 --[GENE_RESPONDS_TO_CONDITION]--> drought
  OsERF83 --[GENE_RESPONDS_TO_CONDITION]--> high salinity
  OsERF83 --[GENE_RESPONDS_TO_CONDITION]--> low temperature
  OsERF83 --[GENE_RESPONDS_TO_CONDITION]--> ABA
- "OsARF1 is expressed in roots" → NOT this family; use GENE_EXPRESSED_IN.

If a statement says a gene is expressed in a tissue under a stress, split only when
both are actually supported:
- condition response: gene responds to the stress;
- static tissue expression: only if baseline/static tissue expression is also stated.

### 4. GENE_EXPRESSED_IN

Static tissue/cell-type expression or subcellular localization, without a treatment
causing the pattern.
The diagnostic question: "Where is this gene/protein normally expressed or localized?"

Include:
- tissue/cell expression: root, leaf, panicle, spikelet, seed, embryo, vascular tissue;
- subcellular localization: nucleus, cytoplasm, chloroplast, membrane, etc.

Examples:
- "OsERF83 is predominantly expressed in vascular tissue" → target="vascular tissue",
  experiment_basis="promoter_GUS; histochemical_staining" if stated.
- "OsERF83-GFP localized to the nucleus" → target="nucleus",
  target_type="subcellular_location", experiment_basis="GFP_fusion_subcellular_localization".
- "OsERF83 showed high expression in roots" → target="root",
  experiment_basis="qRT-PCR; tissue_expression_analysis".

Tissue as experimental context ≠ expression in that tissue.
Put assay context in condition_context or evidence_text instead.

### 5. GENE_ASSOCIATED_WITH_METABOLITE_OR_PATHWAY

Connection to metabolite or metabolic pathway.
Include: starch, sucrose, amylose, flavonoids, lignin, gibberellins, ABA, ROS,
terpenoid pathway, amino acid transport, lignin biosynthesis pathway.

Use this for pathway/metabolite-level claims. If the text lists specific downstream
genes in a pathway, also emit individual GENE_INTERACTS_WITH_GENE_OR_PROTEIN relations.

### 6. GENE_INTERACTS_WITH_GENE_OR_PROTEIN

Physical interaction, genetic interaction, co-expression, transcriptional regulation,
or downstream gene/protein regulation.
Include: protein-protein interaction, TF-target, kinase-substrate, co-expression,
RNA-seq/qRT-PCR downstream expression effects, ChIP/EMSA/luciferase-supported targets.

Downstream gene splitting rule:
When multiple downstream genes are listed by name and the text supports regulation,
extract each named gene as a separate relation. Do not collapse validated named genes
into a single group.

Example:
"OsERF83 upregulated OsNPF8.10, OsNPF8.17 and OsLHT1" → emit three relations.

Groups are allowed only when individual names are unavailable or the text uses a broad
class without listing members, e.g., "PR genes", "WAK genes", "LysM-RLK genes".

For RNA-seq-only regulation, set support_type="evidence" and experiment_basis="RNA-seq".
For RNA-seq plus qRT-PCR validation, set experiment_basis="RNA-seq; qRT-PCR_validation".

### 7. GENE_ASSOCIATED_WITH_VARIANT_QTL_LOCUS

Association with a variant, QTL interval, locus, candidate region, haplotype, or allele.
Include: candidate gene for qGL3, associated with grain width QTL, SNP association,
haplotype-trait association.

### 8. GENE_BELONGS_TO_FAMILY

Use this relation for gene family membership.

Gene family exhaustive expansion rule:
- For gene family papers, extract all explicitly listed or clearly defined family
  members as individual GENE_BELONGS_TO_FAMILY relations.
- If the paper provides a complete numbered range, expand the range into individual
  genes.
  Example: "ONAC001 through ONAC075 were identified as NAC family genes" → emit
  75 individual relations: ONAC001 -> NAC, ONAC002 -> NAC, ..., ONAC075 -> NAC.
- Do not emit compact range relations such as "ONAC001-ONAC075 -> NAC".
- Do not emit GENE_FAMILY_HAS_SIZE.
- Do not emit GENE_FAMILY_IN_SPECIES.
- Family size can be derived later by counting GENE_BELONGS_TO_FAMILY relations.
- Species is stored separately from the caller-provided species field.

If only a few named family members are listed, emit those individual members.
If specific members are experimentally characterized later in the paper, emit those
characterized members individually as well.

Family ID naming: use abbreviation as family_id (ARF, NAC, ERF, bZIP, MADS).
Use the full family name if stated.



Minimal JSON for repeated family membership:
- For routine GENE_BELONGS_TO_FAMILY rows in a large family table, emit ONLY these
  fields if possible: subject_gene, relation_family, target, target_type,
  support_type, evidence_text, experiment_basis.
- Omit optional null/empty fields for these repeated family rows instead of writing
  many null keys. The downstream normalizer accepts missing optional keys.
- Use very short evidence_text, e.g. "Table 1 lists ONAC001-ONAC075 as rice NAC proteins."
- Do not include raw_relation_wording, polarity, condition_context, conclusion_text,
  experiment_methods, analysis_methods, subject_gene_id, target_gene_id for routine
  repeated family rows unless essential.

Compact family-relation fields:
- For repeated GENE_BELONGS_TO_FAMILY relations in a large family table, keep each relation compact.
- Use short shared evidence_text such as "Table 1 lists ONAC001-ONAC075 as rice NAC proteins." Do not repeat long paragraphs.
- Set experiment_basis to a concise term such as "computational_gene_family_identification".
- Do not fill long conclusion_text for routine family membership rows; use null unless needed.
- Do not include gene_id_links for every family member unless the caller specifically asks.

---

## Extraction rules

**Rule 1.** Work only from the provided text. Do not use outside biological knowledge
unless explicitly supported.

**Rule 2.** Prefer explicit assertions. Extract what the document directly states,
not what seems biologically plausible.

**Rule 3.** Preserve raw relation wording in raw_relation_wording as close to the text
as possible: regulates, is induced by, affects, interacts with, encodes, localizes to,
is expressed in, no significant effect, etc.

**Rule 4.** Every relation must have at least one evidence sentence or sentence fragment.

**Rule 5.** Prioritize relations with clear direct evidence and biological relevance.
Do not privilege a relation merely because a gene appears in the title.

**Rule 6.** If no clear evidence exists for any relation, emit no relations for that gene.

**Rule 7.** Distinguish static expression/localization from treatment response:
- GENE_EXPRESSED_IN = where the gene/protein is expressed or localized.
- GENE_RESPONDS_TO_CONDITION = how a condition/treatment changes the gene/protein.

**Rule 8.** Avoid duplicate near-identical relations. But do not under-extract distinct
categories: stress response, tissue expression, localization, phenotype, physiological
readout, and downstream gene regulation are distinct categories.

**Rule 9.** Downstream readout patterns may be extracted when they are direct findings
of this paper. Mark the assay in experiment_basis and keep the evidence precise.

**Rule 9a.** Preserve molecular activity/protein readouts when explicitly linked to a gene
product. If a treatment changes both protein level and associated activity, preserve both
in evidence_text/conclusion_text and emit an activity/function relation when directly supported.

**Rule 9b.** Preserve treatment/condition mechanism clauses in condition_context or
conclusion_text. Do not create a new gene relation solely from a treatment mechanism unless
the text directly connects that mechanism to the gene.

**Rule 9c.** Preserve physiological assay readouts when directly used to support a gene
phenotype. Examples: survival rate, chlorophyll fluorescence, photochemical efficiency,
JIP test, photosynthetic rate, water loss rate, ROS, MDA, proline, electrolyte leakage,
yield components. Use GENE_ASSOCIATED_WITH_TRAIT_OR_PHENOTYPE and include the assay in
experiment_basis or analysis_methods.

**Rule 9d.** Extract negative/no-effect findings when meaningful and experimentally tested.
Use polarity="neutral" and put the no-effect statement in conclusion_text.

**Rule 9e.** Split named downstream genes into separate relations. Keep broad gene groups
only if the paper gives no individual names or intentionally discusses the group as a group.

**Rule 10.** For gene family papers: be comprehensive at the member level.
Expand clearly defined family ranges or complete member lists into individual
GENE_BELONGS_TO_FAMILY relations. Do not output GENE_FAMILY_HAS_SIZE or
GENE_FAMILY_IN_SPECIES.

---

## Gene name ↔ Gene ID mapping rule

Only emit a gene mapping when the text itself makes the correspondence explicit.

Strong enough: OsMFT1 (LOC_Os01g04251), LOC_Os01g04251 designated OsMFT1,
Glyma.05G123400 encodes GmXYZ1, same-sentence explicit apposition.

Do NOT emit a mapping when only the name appears, only the ID appears,
or both appear but the text does not explicitly pair them.
If no explicit mapping is present, return an empty gene_id_links array.

Fields: gene_name, gene_id, mapping_type (prefer explicit_text),
evidence_text, confidence (prefer high only).

---

## support_type

Use ONLY these support_type values:
- evidence: directly supported by an experiment, analysis, or observation in this paper
- background: prior knowledge or cited background, not a finding of this paper
- hypothesis: proposed or inferred without direct establishment
- unclear: evidential status cannot be determined

Do not use custom support_type labels such as "expression", "phenotype", "direct_assay",
"tissue_expression_analysis", or "transcriptome_expression". Put those in experiment_basis,
experiment_methods, or analysis_methods.

---

## experiment_basis guidance

When clearly stated, capture specific methods or evidence basis in experiment_basis,
experiment_methods, analysis_methods, experiment_context, and conclusion_text.

Prefer clear, human-readable terms such as:
- overexpression
- loss_of_function_mutant
- crispr_cas9_editing
- rna_interference
- qRT-PCR
- RNA-seq
- microarray_expression
- qRT-PCR_validation
- promoter_GUS
- histochemical_GUS_staining
- GFP_fusion_subcellular_localization
- subcellular_localization
- tissue_expression_analysis
- drought_phenotyping
- salt_stress_treatment
- cold_stress_treatment
- ABA_treatment
- chlorophyll_fluorescence_JIP_test
- photochemical_efficiency_assay
- survival_rate_assay
- yield_component_analysis
- grain_morphology_analysis
- natural_variant_association
- qtl_or_gwas_association
- protein_interaction_assay
- yeast_two_hybrid
- kinase_activity_assay
- protein_level_assay
- cell_cycle_block
- treatment_timecourse
- gene_coexpression

Use null when no clear basis is stated.

The list above is preferred vocabulary, NOT a closed ontology. If the paper states a
method or assay not listed here, still extract the supported relation and use a concise,
human-readable experiment_basis derived from the paper wording.

Never omit a supported relation only because the method is unfamiliar or absent from
this list. If the relation is supported but the exact method is not stated, set
experiment_basis=null and preserve the available context in evidence_text or
condition_context.

Do not guess or infer an assay name. Avoid phrases such as "likely qRT-PCR" or
"probably northern blot". If the assay is not explicitly stated, use null.
