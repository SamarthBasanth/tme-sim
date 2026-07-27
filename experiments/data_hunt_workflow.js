export const meta = {
  name: 'paired-data-hunt',
  description: 'Wide parallel hunt for OPEN datasets that pair tissue/microenvironment structure with MEASURED drug response or measured spatial drug distribution',
  phases: [
    { title: 'Hunt', detail: '28 agents sweep distinct source families for open paired data' },
    { title: 'Verify', detail: 'confirm URLs resolve, licenses are open, and download what is reachable' },
    { title: 'Rank', detail: 'synthesize a ranked, honest map of what can actually close the validation loop' },
  ],
}

const ROOT = '/Users/samarthbasanth/tme-sim'
const VENV = `${ROOT}/.venv/bin/python`

const FIND_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['angle', 'datasets', 'best_candidate', 'downloaded_anything', 'honest_summary'],
  properties: {
    angle: { type: 'string' },
    datasets: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['name', 'url', 'open_status', 'has_response', 'has_spatial_or_molecular', 'readout', 'why_useful'],
        properties: {
          name: { type: 'string' },
          url: { type: 'string', description: 'exact accession/DOI/download URL that you VERIFIED resolves' },
          open_status: { type: 'string', enum: ['open_direct_download', 'open_registration', 'gated_dua', 'unknown'] },
          has_response: { type: 'string', description: 'what MEASURED drug/treatment response exists (IC50, viability, RANO, survival, imaging), or "none"' },
          has_spatial_or_molecular: { type: 'string', description: 'spatial coordinates? imaging? bulk/sc expression? drug-distribution map?' },
          readout: { type: 'string' },
          why_useful: { type: 'string', description: 'how it could validate a microenvironment drug-transport/response simulator' },
          downloaded_path: { type: 'string', description: 'local path if you downloaded it, else empty' },
        },
      },
    },
    best_candidate: { type: 'string' },
    downloaded_anything: { type: 'boolean' },
    honest_summary: { type: 'string', description: 'blunt assessment: does this angle yield anything that can close the loop?' },
  },
}

// Distinct source families / angles. Each hunter owns one so they do not collide.
const ANGLES = [
  { k: 'maldi_msi_drug', p: 'MALDI / mass-spectrometry IMAGING of DRUG DISTRIBUTION in tumor tissue (measured spatial drug concentration maps). Search METASPACE (metaspace2020.eu, open MSI datasets), PRIDE, MassIVE, Zenodo, GEO, and papers by Tannock/Castellino/Stoeckli/Prideaux on drug MSI in tumors. THIS IS THE HIGHEST-VALUE ANGLE: a measured spatial drug map would let a transport simulator be validated pixel-by-pixel. Look for doxorubicin/TMZ/platinum/erlotinib MSI in tumor sections.' },
  { k: 'msi_repositories', p: 'Open mass-spectrometry-imaging repositories generally: METASPACE public projects, MassIVE, PRIDE, ProteomeXchange, EMBL/Zenodo MSI. Find ANY open MSI dataset containing an administered drug or its metabolite in tumor tissue with spatial coordinates (imzML format).' },
  { k: 'dce_mri_response', p: 'TCIA (The Cancer Imaging Archive) GBM collections with DCE/DSC-MRI PERFUSION imaging AND treatment response/outcome: TCGA-GBM imaging, IvyGAP radiology, UPENN-GBM, RSNA-ASNR-MICCAI BraTS, LUMIERE, QIN-GBM. Perfusion maps are measured per-voxel vascular delivery; paired with RANO/survival this is a real spatial-delivery-vs-outcome dataset. Verify open access + what response labels exist.' },
  { k: 'upenn_lumiere', p: 'UPENN-GBM and LUMIERE longitudinal GBM MRI datasets specifically: are imaging + segmentations + treatment + response/progression openly downloadable? What perfusion sequences? Exact download mechanics (TCIA NBIA, aspera, direct).' },
  { k: 'geo_organoid_screens', p: 'GEO/ArrayExpress for PATIENT-DERIVED GBM ORGANOID or explant drug screens with paired molecular data. Search GEO for glioblastoma organoid + drug + temozolomide, and chase supplementary tables of papers (Jacob 2020, Golebiewska 2020, Loong, Abdullah, Shakya). Need a NUMERIC per-model drug-response table, not figure-level.' },
  { k: 'geo_treated_spatial', p: 'GEO/Zenodo for SPATIAL transcriptomics (Visium/Xenium/MERFISH/GeoMx) of TREATED or RECURRENT (post-temozolomide/post-radiation) glioblastoma, ideally paired vs primary. Treated-vs-untreated spatial tissue with documented therapy is a real in-vivo response signal.' },
  { k: 'geomx_dsp', p: 'NanoString GeoMx DSP / Digital Spatial Profiling datasets of glioma with treatment annotation, on GEO or Zenodo or NanoString public data. Region-of-interest spatial profiling paired with treatment/outcome.' },
  { k: 'imc_codex_treated', p: 'Multiplexed tissue IMAGING (imaging mass cytometry, CODEX, MIBI, Vectra/mIF) of glioma or other tumors WITH treatment response annotation, openly available (Zenodo, IDR - Image Data Resource, Bioimage Archive). Single-cell spatial protein + response.' },
  { k: 'idr_bioimage', p: 'Image Data Resource (idr.openmicroscopy.org) and EBI BioImage Archive: open imaging datasets of tumor tissue with drug treatment and quantified response, especially anything with drug distribution or perivascular gradients.' },
  { k: 'pharmacoscopy_exvivo', p: 'Image-based EX VIVO drug screening on patient tumor samples (pharmacoscopy, Snijder lab; Vladimer; ex vivo functional precision medicine). These have per-cell imaging + measured drug killing. Any open data for glioma or any solid tumor? Check Zenodo/GEO/supplementary.' },
  { k: 'slice_organotypic', p: 'ORGANOTYPIC SLICE CULTURE of human glioma treated with drugs — measured response with preserved tissue architecture (the closest thing to microenvironment + measured kill). Hunt papers + any open data/supplementary quantification.' },
  { k: 'pdx_encyclopedia', p: 'PDX pharmacogenomic resources: Novartis PDX Encyclopedia (PDXE, Gao 2015 supplementary — per-model treatment response + omics), PDXNet, EuroPDX, PDMR (pdmr.cancer.gov open drug-response tables). Does any include GBM/glioma with numeric response + molecular data, openly downloadable?' },
  { k: 'pdmr_hcmi', p: 'NCI PDMR and HCMI specifically: enumerate exactly what is OPEN without an MTA — model metadata, omics, and any drug-response tables. Give exact download URLs and state the gating precisely.' },
  { k: 'depmap_gdsc_ctrp', p: 'Large pharmacogenomic screens with paired omics: DepMap (CRISPR/drug + expression), CTRP, PRISM repurposing, GDSC2, NCI-60, NCI-ALMANAC. Focus on the GLIOMA/GBM subset: how many models, what response metric, and is expression paired and open? These validate potency-vs-molecular, not spatial.' },
  { k: 'hgcc_expression', p: 'HGCC (Human Glioblastoma Cell Culture, hgcc.se) — we already have its drug-response matrix. Now find and download its PAIRED molecular data: GSE72218/GSE72209 expression, any CNV/methylation, and the portal API. Goal: per-culture expression + measured drug response, openly paired — a real molecular->response validation set.' },
  { k: 'cgga_clinical', p: 'CGGA (Chinese Glioma Genome Atlas, cgga.org.cn) — open RNA-seq + clinical with TMZ/radiation treatment and outcome. Also GLASS consortium open tiers, and any glioma cohort with documented per-patient chemotherapy + response/progression. Verify direct downloads.' },
  { k: 'tcga_treatment_response', p: 'TCGA-GBM/LGG treatment + response fields in the GDC: which patients have recorded drug names + treatment_or_therapy + measure_of_response? Quantify how many have a usable RECIST-like or progression response paired with molecular data. Be precise about what response fields actually contain.' },
  { k: 'clinical_trial_open', p: 'Open clinical-trial data-sharing platforms (Vivli, YODA, Project Data Sphere, ClinicalStudyDataRequest) for GLIOMA trials with per-patient treatment + response, and whether any tissue/imaging accompanies it. State the request process precisely.' },
  { k: 'biorxiv_recent', p: 'bioRxiv/medRxiv preprints (2024-2026) announcing OPEN datasets that pair spatial/tissue architecture with measured drug response in ANY cancer (not just GBM) — spatial + perturbation, spatial pharmacology, drug-distribution imaging. Chase their data-availability statements to actual repositories.' },
  { k: 'spatial_perturbation', p: 'Emerging "spatial perturbation / spatial pharmacology" datasets: perturb-seq with spatial readout, spatially-resolved drug response, tumor-on-chip with imaging + kill quantification. Any open data in any cancer that pairs architecture with measured drug effect.' },
  { k: 'organ_on_chip', p: 'Tumor-on-chip / microphysiological-system datasets with measured drug penetration or kill AND geometry (channel/vessel structure). Open data from Emulate/Mimetas/academic MPS groups, Zenodo/GEO. Also microfluidic glioma models.' },
  { k: 'zenodo_figshare', p: 'Zenodo, Figshare, Dryad, OSF: broad sweep for open datasets combining tumor tissue structure/imaging with drug response or drug concentration measurements. Use their search APIs. Report the genuinely useful hits with DOIs.' },
  { k: 'synapse_open', p: 'Synapse.org open-tier projects: GLASS, PS-ON, HTAN open tiers, CPTAC, Brain Tumor projects — what is downloadable WITHOUT a DUA, and does anything carry treatment response paired with tissue/spatial data?' },
  { k: 'htan_open_tier', p: 'HTAN specifically: enumerate exactly which HTAN GBM/brain data is in the OPEN tier (Level 1-4 distinctions, the HTAN Data Portal / CDS), what requires dbGaP+DUA, and whether any open HTAN spatial data carries treatment/response annotation. Precise, not hand-wavy.' },
  { k: 'ebi_arrayexpress', p: 'EBI ArrayExpress/BioStudies, ENA, and the Single Cell Expression Atlas for glioma datasets with drug treatment arms and measured response. Also check the EGA public metadata to identify (gated) candidates worth requesting.' },
  { k: 'sc_drug_perturb', p: 'Single-cell RNA-seq of GBM cells/tumors BEFORE and AFTER drug exposure (temozolomide, targeted agents) — open GEO datasets where the measured outcome is transcriptional/viability change. These give real measured drug effect on real cells with molecular context.' },
  { k: 'radiopharm_pet', p: 'PET imaging of labeled drug/tracer distribution in brain tumors (e.g. 11C-temozolomide, 18F-fluoropyrimidines, radiolabeled agents) — MEASURED in-vivo drug delivery maps in humans. Any open imaging data (TCIA, Zenodo)? This validates delivery physics directly in patients.' },
  { k: 'autoradiography_microdialysis', p: 'Quantitative autoradiography and intratumoral MICRODIALYSIS studies measuring drug concentration inside human/animal brain tumors (e.g. TMZ microdialysis, Portnow-style). Find datasets or precise published measurements (with PMIDs) usable as spatially-resolved or regional PK validation targets.' },
]

function huntPrompt(a) {
  return `You are one hunter in a wide, fast, parallel search for OPEN data that can validate a patient-specific tumor-microenvironment drug simulator at ${ROOT}.

WHAT WE ALREADY HAVE (do NOT re-find these): 10x Xenium GBM single-cell spatial; 4 Visium GBM sections; Darmanis GSE84465 scRNA atlas; Ivy GAP; GDSC glioma; HGCC drug-response matrix; PMID-verified tissue-PK penetration/partition anchors.

WHAT WE CRITICALLY NEED — anything that pairs TISSUE STRUCTURE / MICROENVIRONMENT with a MEASURED DRUG EFFECT. In priority order:
  (1) MEASURED SPATIAL DRUG DISTRIBUTION in tumor tissue (drug concentration vs position) — validates predicted penetration maps directly.
  (2) Tissue/spatial data paired with MEASURED per-sample drug response (kill/IC50/RANO/progression).
  (3) Per-model molecular data paired with measured drug response (no spatial) — validates the molecular->response axis.
  (4) Measured regional/intratumoral drug concentrations (microdialysis, autoradiography, PET) — validates delivery physics.

YOUR ANGLE (own it, don't stray): ${a.p}

METHOD — be fast and rigorous:
- Use WebSearch/WebFetch (load via ToolSearch: "select:WebSearch,WebFetch"). PubMed MCP tools are available too for citations.
- VERIFY every URL you report actually resolves (curl -sI or WebFetch). Do NOT report a link you did not check.
- State open status precisely: open_direct_download / open_registration / gated_dua / unknown.
- If something is OPEN and reachable and reasonably sized (< ~500 MB), DOWNLOAD it into ${ROOT}/data_cache/hunt/${a.k}/ using curl, and record the local path. Verify the file parses/opens (venv python: ${VENV}). Do not download giant archives; prefer the smallest complete unit or a supplementary table.
- Prefer NUMERIC data tables over figure-level claims. A paper's figure is not a dataset.

BE BRUTALLY HONEST. Most angles will come up empty or gated — say so plainly in honest_summary. A precise "nothing open here, and here is exactly what is gated and how to request it" is a valuable result. Do NOT overstate a dataset's pairing.

Return the schema. Report every genuinely relevant dataset you verified (up to ~6), your single best_candidate, and whether you downloaded anything.`
}

phase('Hunt')
const finds = (await parallel(ANGLES.map(a => () =>
  agent(huntPrompt(a), { schema: FIND_SCHEMA, label: `hunt:${a.k}`, phase: 'Hunt' })
))).filter(Boolean)
log(`hunters reported: ${finds.length}/${ANGLES.length}`)

const allDatasets = finds.flatMap(f => (f.datasets || []).map(d => ({ ...d, angle: f.angle })))
const openOnes = allDatasets.filter(d => d.open_status === 'open_direct_download' || d.open_status === 'open_registration')
log(`datasets found: ${allDatasets.length} (open: ${openOnes.length}); downloads attempted by ${finds.filter(f => f.downloaded_anything).length} hunters`)

phase('Verify')
// Independently verify the most promising OPEN candidates actually deliver usable paired data.
const topOpen = openOnes
  .filter(d => (d.has_response || '').toLowerCase() !== 'none')
  .slice(0, 10)
const verified = topOpen.length ? (await parallel(topOpen.map((d, i) => () =>
  agent(`You are independently VERIFYING a candidate dataset for a tumor-microenvironment simulator at ${ROOT}. Be skeptical: hunters overstate pairing.

CANDIDATE: ${d.name}
URL: ${d.url}
Claimed open status: ${d.open_status}
Claimed response data: ${d.has_response}
Claimed spatial/molecular: ${d.has_spatial_or_molecular}
Claimed local download: ${d.downloaded_path || '(none)'}

DO:
1. Actually fetch/inspect it (WebFetch/curl; load tools via ToolSearch "select:WebSearch,WebFetch"). If a local path was claimed, verify the file EXISTS and PARSES with ${VENV}.
2. Answer concretely: does it contain a MEASURED drug/treatment response (numeric, per sample/model/patient)? Does it contain tissue structure (spatial coords, imaging, or drug-concentration-vs-position)? Are those two things linked by a shared identifier — i.e. genuinely PAIRED, not two separate things in one paper?
3. If open and reachable and <500MB, download the smallest usable unit into ${ROOT}/data_cache/hunt/verified_${i}/ and confirm it parses.
4. Report the exact columns/fields that carry the response and the structure.

Return the schema with datasets=[the single verified entry with corrected fields], best_candidate=the name, and an honest_summary that states plainly whether this can close the validation loop (spatial+response), partially close it (molecular+response), or cannot. Do not be generous — say "not genuinely paired" if that is the truth.`,
    { schema: FIND_SCHEMA, label: `verify:${(d.name || 'cand').slice(0, 28)}`, phase: 'Verify' })
))).filter(Boolean) : []
log(`verified ${verified.length} top candidates`)

phase('Rank')
const synthesis = await agent(
  `Synthesize a rigorous, HONEST map of open data that could validate a patient-specific tumor-microenvironment drug-transport/response simulator.

HUNTER RESULTS (${finds.length} angles):
${JSON.stringify(finds, null, 1).slice(0, 60000)}

INDEPENDENT VERIFICATION of top candidates:
${JSON.stringify(verified, null, 1).slice(0, 25000)}

Write ${ROOT}/results/DATA_HUNT.md with:
1. **Headline**: does ANY open dataset pair tissue structure with measured drug response/distribution? Answer plainly.
2. **Tier 1 — can close the loop** (measured spatial drug distribution, OR spatial tissue + measured response): list each with URL, open status, exactly what it contains, and what validation it enables. If empty, say so unambiguously.
3. **Tier 2 — partial** (molecular + measured response, no spatial; or regional/measured PK): these validate the potency or delivery axis. Rank them.
4. **Tier 3 — gated**: what is behind dbGaP/DUA/MTA, and the precise request path (portal, form, expected timeline).
5. **Recommended immediate action**: the single best open dataset to run on NOW, and exactly what analysis it supports.
6. **Honest limits**: what still cannot be validated with open data.
Use only verified URLs. Do not overstate pairing.

Also return a <=350-word summary with the Tier-1 verdict and the top-3 actionable datasets (name + what it enables + local path if downloaded).`,
  { label: 'synthesis', phase: 'Rank' })

return {
  angles: finds.length,
  datasets_found: allDatasets.length,
  open_datasets: openOnes.length,
  downloaded_locally: allDatasets.filter(d => d.downloaded_path).map(d => ({ name: d.name, path: d.downloaded_path })),
  verified_count: verified.length,
  synthesis,
}
