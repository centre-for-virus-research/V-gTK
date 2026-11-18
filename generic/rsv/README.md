# Human RSV (*Orthopneumovirus hominis*) Reference Set (`generic/rsv/`)

- **NCBI Taxonomy ID**: `3049954` (*Orthopneumovirus hominis*, species rank; encompasses `11250` Human RSV, `208893` RSV-A, `208895` RSV-B, and modern species-level submissions).
- **Genome architecture**: Non-segmented (`is_segmented = "N"`), negative-sense ssRNA (~15.2 kb), 11 canonical ORFs (`NS1`, `NS2`, `N`, `P`, `M`, `SH`, `G`, `F`, `M2-1`, `M2-2`, `L`).

## Reference Curation (`ref_list.txt` — 30 complete genomes)

All 30 accessions are bare primary INSDC/RefSeq accessions verified inside `txid3049954[Organism:exp]`, complete genomes (`15,111–15,283 nt`), and free of `N` stretches (`28/30` have zero ambiguous bases; `2/30` have a single IUPAC base):

1. **Master reference (`1` row)**:
   - `NC_038235` (`A / GA1`, `15,222 nt`): NCBI RefSeq for Human RSV-A with all **11 CDSs** annotated on `+` strand (unlike older `NC_001803`, which omits `M2-2`). `F` spans `5661..7385` (`574` codons + stop).
2. **Subgroup A references (`15` rows)**:
   - `GA1` classic prototypes (`AY911262` = ATCC VR-26 **Long** strain, focal background of Simonich et al. 2026 *Nature* `RSV_Long_F_DMS`; `KT992094` & `M74568` = **A2** structural prototype; `MG642074` = 1980 clinical `GA1`).
   - Pre-duplication clinical lineages `GA5` (`JX015498`, `KF826849`) and `NA1` (`JF920069`, `KJ627352`).
   - Modern `ON1` / `A.D.*` lineages carrying the **72-nt `G`-gene duplication** (`KC731482`, `JX627336`, `KY654518`, `OK649682`, `ON729320`, `OP890319`, `OZ556678` — 7 diverse 2011–2025 genomes, well above `ref_aln_min_insertion_support = 2` so the 72-nt `G` duplication is preserved in-frame in the reference backbone).
3. **Subgroup B references (`14` rows)**:
   - Non-duplicated `GB1`–`GB3` prototypes (`NC_001781` = NCBI RSV-B RefSeq B1; `AF013254` = wildtype B1; `AY353550` = strain 9320; `KP258736` = 1982 `GB1`; `MG642066` = `GB2`; `KF826829` = 2005 `GB3`).
   - Modern `BA` / `BA9` / `B.D.*` lineages carrying the **60-nt `G`-gene duplication** (`JQ582844`, `KF640637`, `KX765890`, `KY249660`, `OP975389` [Nextstrain RSV-B reference `VIC-RCH056/2019`], `OP965709`, `LC817403`, `OZ556688` — 8 diverse 2004–2025 genomes so the 60-nt `G` duplication is preserved in-frame in the reference backbone).

