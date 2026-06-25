# Literature Validation: High Neutrophil Proportions in Captive Bat Whole Blood

**Context:** scRNA-seq of captive *Eonycteris spelaea* (cave nectar bat) whole blood. After correcting a cross-species annotation artifact, neutrophils = 0% in some samples but 25–63% in others (ES18 63%, ES171 59%). Question: are high neutrophil fractions biologically plausible in bat blood?

**Bottom line up front:** 25–63% neutrophils in captive bat whole blood is biologically credible and well within / near published pteropodid and general mammalian ranges. The 0% values are the artifact (annotation + known scRNA-seq granulocyte dropout), not the high values. NO direct *E. spelaea* peripheral-blood neutrophil-% reference range was found in the literature; the closest data come from other pteropodid/fruit bats and from the *E. spelaea* and Egyptian rousette single-cell studies (which report neutrophils as a normal, sometimes dominant, myeloid population).

---

## Point 1 — Eonycteris spelaea / Pteropodid (fruit/nectar bat) blood leukocyte differential

**No published peripheral-blood neutrophil reference *percentage* specific to *E. spelaea* was located.** The flagship single-cell study does not report a clinical differential, but it explicitly identifies and characterizes a normal bat neutrophil population:

- **Gamage AM, Chan WOY, Zhu F, … Wang LF. (2022).** "Single-cell transcriptome analysis of the in vivo response to viral infection in the cave nectar bat *Eonycteris spelaea*." *Immunity* 55(11):2187–2205.e5. DOI: 10.1016/j.immuni.2022.10.008. PMID 36351376.
  - This IS the *E. spelaea* single-cell immune atlas (Aaron Irving / Lin-Fa Wang group, Duke-NUS). It defines bat neutrophils as a distinct population ("Bat neutrophils were distinguished by high basal IDO1 expression"). Confirms neutrophils are a genuine, expected immune compartment in *E. spelaea*; their presence is normal, not anomalous. (Tissue analyzed was lung, not a clinical blood differential — so no % blood reference is given here.)
  - Companion commentary: **Banerjee A, Mossman K. (2022).** "Laying the foundation for single-cell studies in bats." *Immunity* 55(11):1974–1977. DOI: 10.1016/j.immuni.2022.10.010. PMID 36351371.

- **Friedrichs V, Toussaint C, Schäfer A, … Saliba AE, Dorhoi A. (2022).** "Landscape and age dynamics of immune cells in the Egyptian rousette bat." *Cell Reports* 40(10):111305. DOI: 10.1016/j.celrep.2022.111305. PMID 36070695.
  - scRNA-seq + flow cytometry of *Rousettus aegyptiacus* (pteropodid fruit bat). Key finding: **"neutrophils, CD206+ myeloid cells, and CD3+ T cells dominate as bats reach adulthood"** — i.e., neutrophils are a DOMINANT leukocyte population in adult fruit bats, with a predominance of neutrophils in adults vs. lymphocytes in juveniles. Directly supports high neutrophil fractions in adult pteropodids. (Cross-references *E. spelaea* and *Pteropus* in its comparative discussion.)

**Pteropodid / fruit-bat hematology reference studies (closest available numeric data):**
- **Selig M, Lewandowski A, Kent MS. (2016).** "Establishment of reference intervals for hematology and biochemistry analytes in a captive colony of straw-colored fruit bats (*Eidolon helvum*)." *J Zoo Wildl Med* 47(1):106–112. DOI: 10.1638/2015-0040.1. PMID 27010270. (n=45 captive megachiropterans; values stated to be similar to other pteropodid bat species.)
- **Ekeolu OK, Adebiyi OE. (2018).** "Hematology and erythrocyte osmotic fragility of the Franquet's fruit bat (*Epomops franqueti*)." *J Basic Clin Physiol Pharmacol* 29(4):391–394. DOI: 10.1515/jbcpp-2017-0169. PMID 29543589. (n=17; full differential WBC including neutrophil count reported; reference-range establishment study.)
- **Strumpf AA, Malmlov A, Ayers JD, Schountz T, Kendall LV. (2020).** "Hematologic values of Jamaican fruit bats (*Artibeus jamaicensis*) and the effects of isoflurane anesthesia." *J Am Assoc Lab Anim Sci* 59(3):275–281. DOI: 10.30802/AALAS-JAALAS-19-000056. PMID 32164795; PMCID PMC7210728. (Captive colony baseline differential; also a handling/anesthesia effect study — see Point 4.)

> NOTE: Bat granulocytes are frequently termed "neutrophils" in megachiropterans (true neutrophils) but "heterophils" is used in some chiropteran/avian-style reports. Both denote the same granulocyte compartment for proportion purposes.

---

## Point 2 — Mammalian whole-blood neutrophil fraction (context that >50% is ordinary)

In human whole blood, neutrophils are the most abundant leukocyte, normally **~40–70% (commonly ~50–70%) of circulating leukocytes** — standard clinical hematology (absolute ~1.8–7.7 ×10⁹/L). Many domestic mammals (dog, horse, primates) are similarly neutrophil-predominant. Therefore a 50–63% neutrophil fraction in whole blood is the textbook-normal state for a neutrophil-dominant mammal, not an outlier.

- Comparative immune-proportion framing across taxa: **Cornelius Ruhs E, Becker DJ, Oakey SJ, … Downs CJ. (2021).** "Body size affects immune cell proportions in birds and non-volant mammals, but not bats." *J Exp Biol* 224(13):jeb241109. DOI: 10.1242/jeb.241109. PMID 34104965. Quantifies lymphocyte vs granulocyte (neutrophil/heterophil) proportions across **63 bat species, 400 bird species, 251 non-volant mammal species** — establishes that granulocyte (neutrophil) proportions are a major, variable fraction of leukocytes across bats and mammals; bat granulocyte proportions overlap the mammalian range. (Confirms neutrophils/granulocytes form a large, normal share of the leukogram in bats.)

---

## Point 3 — Neutrophil/granulocyte dropout in droplet scRNA-seq (so 0% is the artifact)

Strong, well-established evidence that granulocytes (neutrophils, eosinophils) are systematically UNDER-captured or entirely ABSENT from droplet/microwell scRNA-seq atlases due to low mRNA content, high RNase, and fragility — meaning scRNA-seq proportions UNDERSTATE the true neutrophil fraction and can drop it to ~0:

- **Borrelli C, Gurtner A, Arnold IC, Moor AE. (2024).** "Stress-free single-cell transcriptomic profiling and functional genomics of murine eosinophils." *Nat Protoc* 19(6):1679–1709. DOI: 10.1038/s41596-024-00967-3. PMID 38504138.
  - Explicit: granulocytes are **"absent from conventional single-cell RNA sequencing atlases owing to technical difficulties preventing their transcriptomic interrogation"**; requires special low-shear, fast, RNase-minimizing protocols. States the approach "can be adapted to investigate other granulocytes, such as **neutrophils**." Direct support that 0% granulocytes is a known technical artifact.

- **Grieshaber-Bouyer R, Radtke FA, Cunin P, … Nigrovic PA; ImmGen Consortium. (2021).** "The neutrotime transcriptional signature defines a single continuum of neutrophils across biological compartments." *Nat Commun* 12(1):2856. DOI: 10.1038/s41467-021-22973-9. PMID 34001893.
  - Authoritative ImmGen neutrophil single-cell reference; documents the special handling/low-RNA challenges of profiling neutrophils by scRNA-seq (neutrophils are low-RNA, fragile, easily lost in standard pipelines).

- **Vafadarnejad E, Rizzo G, Krampert L, … Saliba AE, Cochain C. (2020).** "Dynamics of cardiac neutrophil diversity in murine myocardial infarction." *Circ Res* 127(9):e232–e249. DOI: 10.1161/CIRCRESAHA.120.317200. PMID 32811295.
  - Demonstrates neutrophil scRNA-seq requires dedicated enrichment/handling; reinforces that standard droplet workflows lose neutrophils.

**Interpretation:** In a whole-blood 10x-style run, observing 0% neutrophils in some samples is the expected failure mode of granulocyte dropout; samples that DID retain 25–63% neutrophils more faithfully reflect the true (high) blood neutrophil fraction.

---

## Point 4 — Stress neutrophilia / handling leukogram (relevant: captive bats, handled at blood draw)

Acute handling/restraint and capture stress reliably cause neutrophilia / elevated neutrophil(heterophil):lymphocyte (N:L or H:L) ratio via catecholamine-driven demargination and glucocorticoid effects, in birds, mammals, and bats:

- **Strumpf AA, et al. (2020).** *J Am Assoc Lab Anim Sci* 59(3):275–281. PMID 32164795 (full cite above).
  - In Jamaican fruit bats, the method of handling (physical restraint vs isoflurane) significantly changed the leukogram — **significant increases in WBC, lymphocytes, and monocytes** between handling conditions — directly demonstrating that handling/restraint perturbs the bat leukocyte differential. Establishes that blood-draw handling is a real confounder of bat differentials.

- **Pusch EA, Bentz AB, Becker DJ, Navara KJ. (2018).** "Behavioral phenotype predicts physiological responses to chronic stress in proactive and reactive birds." *Gen Comp Endocrinol* 255:71–77. DOI: 10.1016/j.ygcen.2017.10.008. PMID 29051076.
  - Shows **heterophil:lymphocyte (H:L) ratios significantly elevated by stress** alongside corticosterone — the canonical stress-leukogram pattern (granulocyte rise relative to lymphocytes). H:L ratio is the standard avian/reptile/wildlife stress biomarker; the chiropteran analogue is N:L ratio.

- **Becker DJ, Dyer KE, Lock LR, Fenton MB, Simmons NB. (2025).** "Habitat and seasonal drivers of leukocyte profiles within and across Neotropical bat species." *Biol Lett* 21(10):20250447. DOI: 10.1098/rsbl.2025.0447. PMID 41159440; PMCID PMC12570079.
  - Uses **neutrophil-to-lymphocyte ratio (N:L)** as a haematological stress/inflammation marker in wild bats; documents strong species-specific variation in total leukocyte counts and N:L. Confirms N:L is the recognized bat stress index and that neutrophil proportions vary widely and respond to physiological state.

N:L ratio as a generalized vertebrate stress marker is further supported in other wild mammals (e.g., **Iwińska K, et al. 2024,** *Biol Lett* 20(10):20240257, DOI 10.1098/rsbl.2024.0257, PMID 39471836; **Beaman JE, et al. 2023,** *Conserv Physiol* 11(1):coac088, DOI 10.1093/conphys/coac088 — koala translocation stress via N:L).

---

## Point 5 — Captive vs wild bat immune/hematology differences

- **Cornelius Ruhs E, et al. (2021).** *J Exp Biol* 224:jeb241109. PMID 34104965 (full cite above).
  - Directly compares **wild-caught vs captive bats** (and captive non-volant mammals/birds) for leukocyte proportions; shows captivity status and ecological context shape granulocyte/lymphocyte proportions, with bats showing distinctive (non-body-size-dependent) granulocyte scaling. Establishes captivity as a modifier of the leukogram.

- **Becker DJ, et al. (2025).** *Biol Lett* 21:20250447. PMID 41159440 (above).
  - Habitat/environmental conditions drive species-specific shifts in total leukocyte counts and N:L ratio in bats — by extension, captive husbandry conditions are expected to alter leukocyte composition and N:L.

- **Strumpf AA, et al. (2020).** PMID 32164795 — provides captive-colony baseline differentials and shows handling sensitivity, the practical captive-bat reference.

**Captivity-since-birth + handling at blood draw** would both be expected to push the leukogram toward higher neutrophil/N:L (chronic captivity stress + acute handling neutrophilia), supporting the observed high neutrophil fractions.

---

## Verdict

**25–63% neutrophils in captive bat whole blood is biologically credible.** It sits within the normal neutrophil-predominant mammalian range (~40–70%, Point 2), matches the finding that neutrophils DOMINATE the adult pteropodid leukogram (Friedrichs 2022; Gamage 2022 confirms a normal *E. spelaea* neutrophil population), and is further expected to be elevated by captivity + handling stress (Points 4–5). Conversely, the 0% samples are explained by the well-documented systematic dropout of granulocytes in droplet scRNA-seq (Point 3) — the 0% is the artifact, the 25–63% is closer to truth. **Caveat:** no *E. spelaea*-specific blood neutrophil reference percentage exists in the literature; the strongest direct support is the Egyptian rousette (Friedrichs 2022, "neutrophils dominate adults") and the general pteropodid + mammalian ranges.
