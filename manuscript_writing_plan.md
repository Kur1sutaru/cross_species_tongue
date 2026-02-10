# Manuscript Writing Plan: Cross-Species Single-Cell Atlas of Tongue
## Updated with New Conservation Analysis Results

---

## CURRENT STATUS ASSESSMENT

### What You Have (Existing Document):
✅ **Strong foundation in place:**
- Title and abstract highlights framework
- Introduction/background with good motivation
- Complete Methods section (ethical approval, data acquisition, computational pipeline)
- Results section with species-specific analyses for each cell type
- Curated marker panels for endothelial, epithelial, fibroblasts, immune, muscle, Schwann cells
- Extensive reference list (21+ citations)

### What's NEW from Your Presentation:
🆕 **Critical updates needed:**
- **Dataset expansion**: Now 135,736 cells (was 76,210)
- **Conservation analysis**: 4,951 conserved genes across 25/32 cell types
- **Statistical framework**: FindConservedMarkers with meta-analysis approach
- **Jaccard index analysis**: Quantified conservation levels (1-3%)
- **Methodological justification**: Why your approach vs alternatives
- **7 major cell types confirmed, 32 subtypes total**
- **Species distribution**: Human 44%, Mouse 35%, Marmoset 15%, Rat 6%

---

## MANUSCRIPT WRITING PLAN

### PHASE 1: UPDATE RESULTS SECTION (Priority: HIGH)
**Timeline: Week 1**

#### 1.1 Create New "Conservation Analysis" Results Section
**INSERT AFTER**: Current section 4.2 (Integration results)
**NEW SECTION**: 4.3 Cross-Species Gene Conservation Analysis

**Content to write:**
```
4.3 Cross-Species Gene Conservation Reveals Core Mammalian Programs

To systematically identify genes that define conserved cell-type identity across 
evolutionary time scales, we performed a rigorous meta-analysis of differential 
expression across all four species.

4.3.1 Conservation Analysis Strategy
[Explain the FindConservedMarkers approach]
- Why independent testing per species (avoid bias toward abundant species)
- Statistical framework (Wilcoxon rank-sum test, meta-p-value combination)
- Thresholds and criteria

4.3.2 Overview of Conserved Gene Landscape
- Successfully analyzed 25/32 cell types (78%)
- 7 excluded (absence in some species or low cell counts)
- Total: 4,951 conserved genes identified
- Average: 190 genes per cell type (range: 56-440)
- Distribution across cell types [create new figure]

4.3.3 Interpreting Conservation Metrics
[Critical: Frame the Jaccard indices positively]
- Jaccard scores 1-3% represent 1,000-3,000x enrichment over random
- Low percentages reflect evolutionary distance, not poor conservation
- These are the genes that truly define cell identity

4.3.4 Cell-Type Specific Conservation Patterns
[For top 5-7 cell types, describe:]
- Number of conserved genes
- Key conserved pathways
- Examples of highly conserved markers
- Biological interpretation
```

**NEW FIGURES NEEDED:**
- **Figure 3A**: Bar plot showing number of conserved genes per cell type (25 cell types)
- **Figure 3B**: Heatmap of top conserved genes across species for select cell types
- **Figure 3C**: Jaccard index analysis with statistical enrichment
- **Figure 3D**: Venn diagrams showing overlap for 2-3 example cell types

#### 1.2 Reorganize Existing Results Sections
**Current structure** → **New structure**

```
OLD:
4.1 Quality control and integration
4.2 Cross-species integration reveals conserved architecture
4.3 Species-specific transcriptional features
4.4-4.9 Individual cell type analyses (endothelial, epithelial, etc.)

NEW:
4.1 Single-cell profiling across four mammalian species
    - Dataset overview (135,736 cells, 7 major types, 32 subtypes)
    - Quality control metrics
    - Species distribution

4.2 Cross-species integration reveals conserved cellular architecture
    - Harmony integration results
    - UMAP visualization
    - Six major cell populations identified

4.3 Cross-species gene conservation analysis [NEW SECTION]
    - Conservation analysis strategy
    - Overview of 4,951 conserved genes
    - Interpreting conservation metrics
    - Cell-type specific patterns

4.4 Species-specific transcriptional divergence
    - Integrate existing cell-type specific lollipop analyses
    - Frame as complementary to conservation analysis
    - Emphasize: conservation + divergence = complete picture
```

#### 1.3 Update Dataset Numbers Throughout
**Global find-replace needed:**
- Total cells: 76,210 → 135,736
- Human cells: [old number] → 59,526 (44%)
- Mouse cells: 47,823 → 47,823 (35%)
- Marmoset cells: 20,063 → 20,063 (15%)
- Rat cells: 8,324 → 8,324 (6%)

---

### PHASE 2: UPDATE METHODS SECTION (Priority: HIGH)
**Timeline: Week 1-2**

#### 2.1 Add Conservation Analysis Methods
**NEW SUBSECTION**: 3.7 Cross-Species Conservation Analysis

**Content to write:**
```
3.7 Cross-Species Conservation Analysis

To identify genes with conserved cell-type-specific expression across species, 
we employed the FindConservedMarkers function in Seurat v5.3.0. This approach 
performs independent differential expression testing within each species and 
combines the results using meta-analysis.

For each cell type:
1. Required presence in all 4 species (minimum 10 cells per species)
2. Identified genes differentially expressed vs. all other cells
3. Statistical testing: Wilcoxon rank-sum test per species independently
4. Meta-p-value: Combined across species using minimump method
5. Thresholds: p-value < 0.05 (Bonferroni-corrected), log2FC ≥ 0.25, 
   minimum detection 25% of cells per species

Rationale for this approach:
- Avoids bias toward abundant species (vs. pooling all species)
- Ensures both high expression AND cross-species consistency (vs. correlation alone)
- Provides statistical rigor with proper ortholog mapping (vs. simple overlap)

The conservative nature of this method (requiring presence in ALL species) 
prioritizes specificity over sensitivity, identifying the most robust and 
evolutionarily conserved cell-type markers.

Jaccard indices were calculated to quantify overlap between species-specific 
marker sets, and statistical enrichment was assessed using hypergeometric tests 
relative to random gene set expectations.
```

#### 2.2 Update Computational Analysis Section
Add mention of conservation analysis tools and parameters to section 3.4

---

### PHASE 3: REVISE ABSTRACT (Priority: MEDIUM)
**Timeline: Week 2**

#### Current Abstract Highlights → New Abstract Structure

**NEW ABSTRACT** (250 words max):

```
[BACKGROUND - 2-3 sentences]
The tongue is a specialized organ comprising diverse cell types essential for 
feeding, speech, and sensory perception. Understanding the evolutionary 
conservation and divergence of tongue cellular architecture is critical for 
translating findings across model systems and identifying therapeutic targets.

[OBJECTIVE - 1 sentence]
We generated a comprehensive cross-species single-cell RNA-seq atlas to define 
conserved and species-specific cell types and transcriptional programs.

[METHODS - 2-3 sentences]
We profiled 135,736 single cells from tongue tissue across four mammalian 
species (human, mouse, marmoset, rat), identifying 32 distinct cell subtypes 
within 7 major lineages. Using ortholog mapping and cross-species integration, 
we performed meta-analysis to identify conserved cell-type markers across all 
species.

[RESULTS - 4-5 sentences]
Cross-species integration revealed highly conserved cellular architecture, 
with epithelial cells, immune cells, fibroblasts, endothelial cells, muscle 
cells, and neural populations present across all species. Conservation analysis 
identified 4,951 genes with statistically significant enrichment in 25 cell 
types across all four species (p < 0.0001), representing core mammalian 
programs maintained over 100+ million years of evolution. Cell-type-specific 
analyses revealed both deeply conserved transcriptional modules (e.g., 
myelinating Schwann cells, smooth muscle) and species-specific adaptations 
(e.g., cornification programs, capillary metabolism, immune activation states).

[CONCLUSIONS - 2 sentences]
This cross-species atlas provides a foundational resource for understanding 
tongue biology, facilitating cross-species translation in oral disease research 
and cancer-neuron interaction studies.

[KEYWORDS]
Single-cell RNA-seq, comparative genomics, tongue, evolutionary conservation, 
cross-species integration, oral cavity, mammalian evolution
```

---

### PHASE 4: STRENGTHEN INTRODUCTION (Priority: MEDIUM)
**Timeline: Week 2**

#### Add Conservation Analysis Motivation

**INSERT PARAGRAPH** (after existing intro paragraphs, before study objectives):

```
While species-specific characterizations provide valuable insights, systematic 
cross-species comparison is essential to distinguish evolutionarily conserved 
core programs from species-specific adaptations. Conserved genes and pathways 
represent fundamental biological principles that are likely translatable across 
model systems, whereas divergent features may reflect species-specific 
adaptations to distinct ecological niches, anatomical constraints, or 
physiological demands [cite relevant cross-species scRNA-seq papers]. Previous 
cross-species single-cell studies in other organs have revealed that while 
cellular identities are broadly conserved, the specific genes defining those 
identities can vary substantially [cite examples]. Therefore, rigorous 
statistical approaches that account for evolutionary distance while identifying 
truly conserved programs are critical for meaningful cross-species inference.
```

**UPDATE STUDY OBJECTIVES** to include:

```
Therefore, in the current study we:
1. Generated single-cell profiles from male and female tongue tissue across 
   four mammalian species (human, mouse, marmoset, rat)
2. Performed cross-species integration to map conserved cellular architecture
3. Identified core conserved genes defining cell-type identity across species [NEW]
4. Characterized species-specific transcriptional divergence
5. Established a resource for cross-species translation in oral biology
```

---

### PHASE 5: EXPAND DISCUSSION (Priority: MEDIUM)
**Timeline: Week 2-3**

#### 5.1 Add Conservation Analysis Discussion

**NEW DISCUSSION SECTIONS:**

```
5.1 Evolutionary Conservation of Tongue Cellular Programs

Our conservation analysis revealed that despite 100+ million years of mammalian 
evolution, tongue cell types retain remarkably stable transcriptional identities. 
The identification of 4,951 conserved genes across 25 cell types demonstrates 
that fundamental cellular programs are preserved even as species diverge in 
morphology, physiology, and behavior.

[Discuss:]
- Interpretation of conservation levels (why 1-3% Jaccard is meaningful)
- Comparison to other cross-species single-cell studies
- Which cell types show highest conservation and why
- Biological implications of conserved programs

5.2 Species-Specific Adaptations Reflect Functional Specialization

While core programs are conserved, we identified substantial species-specific 
variation in expression magnitude and pathway enrichment...

[Integrate existing species-specific findings into evolutionary framework]

5.3 Methodological Considerations for Cross-Species Analysis

Our approach prioritized statistical rigor and specificity over sensitivity...

[Discuss strengths and limitations of FindConservedMarkers approach]

5.4 Implications for Cross-Species Translation

The clear delineation of conserved vs. divergent programs has direct 
implications for translating findings between model systems...

[Discuss translational relevance, which findings are likely conserved vs. 
species-specific]

5.5 Future Directions

[Include spatial transcriptomics, functional validation, disease contexts, 
sex differences, etc.]
```

---

### PHASE 6: CREATE NEW FIGURES (Priority: HIGH)
**Timeline: Week 1-3**

#### Figure Plan

**Figure 1: Study Overview and Dataset** (Update existing)
- Panel A: Experimental workflow (keep)
- Panel B: **UPDATE** - Cell count distribution (135,736 cells, 4 species)
- Panel C: **NEW** - Overview of 7 major cell types, 32 subtypes
- Panel D: QC metrics per species

**Figure 2: Cross-Species Integration** (Keep with minor updates)
- Panel A: UMAP colored by species
- Panel B: UMAP colored by cell type
- Panel C: Integration quality metrics

**Figure 3: Conservation Analysis [NEW FIGURE]** ⭐
- Panel A: Bar plot - # conserved genes per cell type (25 cell types)
- Panel B: Heatmap - top 50 conserved genes across species for 5-6 example cell types
- Panel C: Jaccard index distribution with statistical enrichment
- Panel D: Example Venn diagrams for 3 cell types showing gene overlap

**Figure 4: Cell-Type Specific Conservation and Divergence [NEW INTEGRATED FIGURE]**
- Combine your existing lollipop plots into multi-panel figure
- Panel A: Endothelial markers
- Panel B: Epithelial markers
- Panel C: Fibroblast markers
- Panel D: Immune markers
- Panel E: Muscle markers
- Panel F: Schwann cell markers

**Figure 5-X: Additional cell-type specific analyses**
- Keep existing detailed analyses as supplementary or individual figures

**Supplementary Figures:**
- Individual cell type UMAPs
- QC plots
- Additional marker analyses
- Pathway enrichment for each cell type

---

### PHASE 7: SUPPLEMENTARY MATERIALS (Priority: LOW)
**Timeline: Week 3-4**

#### Supplementary Tables Needed:

**Table S1:** Dataset statistics per species (cells, genes, QC metrics)

**Table S2:** Complete list of 4,951 conserved genes with:
- Gene symbol
- Cell type
- Species-specific log2FC
- Meta-p-value
- Average expression per species

**Table S3:** Cell type marker genes (all 32 subtypes)

**Table S4:** Conservation analysis parameters and thresholds

**Table S5:** Pathway enrichment results for conserved genes

**Table S6:** Species-specific marker genes

---

## WRITING WORKFLOW RECOMMENDATIONS

### Week 1: Core Results & Methods
**Priority tasks:**
1. ✍️ Write new conservation analysis results section (4.3)
2. ✍️ Add conservation methods subsection (3.7)
3. 📊 Create Figure 3 (conservation analysis)
4. 🔄 Update all dataset numbers throughout manuscript
5. 📊 Generate Supplementary Table S2 (conserved gene list)

### Week 2: Structure & Integration
**Priority tasks:**
6. ✍️ Revise abstract with new findings
7. ✍️ Expand introduction with conservation motivation
8. 🔄 Reorganize results section structure
9. ✍️ Start discussion section on conservation
10. 📊 Update Figure 1 with new dataset overview

### Week 3: Discussion & Figures
**Priority tasks:**
11. ✍️ Complete discussion section (all subsections)
12. 📊 Create Figure 4 (integrated lollipop plots)
13. ✍️ Write figure legends for all new figures
14. 🔄 Integrate existing species-specific analyses into new framework
15. 📊 Finalize supplementary figures

### Week 4: Polish & Supplementary
**Priority tasks:**
16. ✍️ Complete all supplementary tables
17. 📝 Proofread entire manuscript for consistency
18. 🔍 Check all cross-references (figures, tables, citations)
19. 📊 Format all figures to journal specifications
20. 📋 Prepare cover letter and highlights

---

## KEY MESSAGES TO EMPHASIZE

### Scientific Contributions:
1. **Scale**: Largest cross-species tongue atlas to date (135,736 cells, 4 species)
2. **Rigor**: Statistically robust conservation analysis identifying 4,951 core genes
3. **Completeness**: Both conservation AND divergence characterized
4. **Resource**: Foundational atlas for oral biology and cancer-neuron research

### Methodological Advances:
1. Rigorous approach avoiding common pitfalls (species bias, correlation vs. expression)
2. Proper statistical framework with meta-analysis
3. Clear interpretation of conservation metrics
4. Comprehensive cell-type coverage (32 subtypes)

### Biological Insights:
1. Core mammalian programs maintained over evolutionary time
2. Species-specific adaptations reflect functional specialization
3. Translational implications clearly delineated
4. Foundation for disease studies and therapeutic target discovery

---

## INTEGRATION STRATEGY

### How to Merge Old + New Content:

#### ✅ KEEP (from existing document):
- All methods sections (with additions)
- Species-specific lollipop analyses (reframe as divergence analysis)
- Introduction biological background
- Reference list (add conservation papers)

#### 🆕 ADD (from presentation):
- Conservation analysis results (new section 4.3)
- Conservation methods (new section 3.7)
- Updated dataset statistics throughout
- New Figure 3 (conservation analysis)
- Discussion of conservation implications

#### 🔄 REVISE (combine old + new):
- Abstract (incorporate conservation findings)
- Introduction (add conservation motivation)
- Results structure (reorganize for better flow)
- Discussion (frame species-specific as complementary to conservation)

---

## NEXT STEPS CHECKLIST

### Before You Start Writing:
- [ ] Finalize exact numbers for all dataset statistics
- [ ] Generate complete conserved gene list (Table S2)
- [ ] Create Figure 3 panels (conservation analysis)
- [ ] Review 3-5 key cross-species single-cell papers for citation
- [ ] Decide on target journal (impacts format, length, style)

### As You Write:
- [ ] Maintain consistent terminology (conserved genes, not preserved/retained)
- [ ] Use consistent species abbreviations (H. sapiens, M. musculus, etc.)
- [ ] Number all new figures and tables consecutively
- [ ] Add internal cross-references as you write
- [ ] Keep a running list of citations needed

### Quality Control:
- [ ] Check that all claims are supported by data/figures
- [ ] Verify all statistics are correctly reported
- [ ] Ensure figure legends are complete and self-explanatory
- [ ] Confirm all methods are adequately described
- [ ] Review for logical flow and coherence

---

## POTENTIAL JOURNAL TARGETS

Based on scope and impact:

**Tier 1 (High Impact):**
- Nature Communications (broad audience, cross-species emphasis)
- Cell Reports (cell biology focus)
- Genome Biology (methods + resource)

**Tier 2 (Specialized):**
- PLOS Biology (open access, broad scope)
- eLife (evolutionary + systems biology)
- Nucleic Acids Research (database + resource)

**Tier 3 (Focused):**
- BMC Biology
- Frontiers in Cell and Developmental Biology
- Scientific Data (if framed as resource paper)

---

## FINAL NOTES

### Strengths of Your Study:
✅ Largest sample size in the field
✅ Rigorous statistical approach
✅ Both conservation and divergence characterized
✅ Multiple species (not just mouse vs human)
✅ Complete methods already documented
✅ Sex differences data available
✅ Clear translational relevance

### What Makes This Manuscript Strong:
The integration of conservation analysis with species-specific characterization 
provides a complete picture. Most papers do one or the other - you're doing both.

The statistical rigor (meta-analysis, proper ortholog mapping, conservative 
thresholds) will satisfy reviewers.

The scale (135K cells, 4 species, 32 subtypes) is impressive and represents a 
significant resource for the field.

### Timeline Estimate:
- **3-4 weeks**: Complete first draft with all new content
- **1-2 weeks**: Internal review and revision
- **1 week**: Finalize figures and formatting
- **Total: 5-7 weeks to submission**

Good luck with the manuscript! The data is strong and the story is clear. 
You have all the pieces - now it's about assembling them into a compelling narrative.
