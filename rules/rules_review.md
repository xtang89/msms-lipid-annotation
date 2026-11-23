# Rule notes (v0.7.0, 2026-11-22)

The rule files (`positive.yaml`, `negative.yaml`) hold only general rules: class, adduct, fragment formula and role.
Each rule is confirmed in the literature (`references.md`). Rules that we could not confirm are not included (table below).

## Scoring
Each fragment has a role:
- **diagnostic**: specific to the head group with this adduct; one is enough to identify the class.
- **supporting**: agrees with the class; increases confidence; cannot identify the class alone.
- **chain**: gives the fatty-acid chains; for TG, DG and MG (no head group) it is also the class evidence.

| Level | Condition | Result |
|---|---|---|
| high | >= 2 diagnostic, or 1 diagnostic + >= 1 supporting | MS/MS confirmed |
| medium | 1 diagnostic only | MS/MS confirmed |
| low | supporting only | MS1 annotation; MS/MS does not confirm the class |
| ms1_only | no matched fragment | MS1 annotation |

- TG, DG, MG: >= 2 chain fragments that fit one chain combination with the correct sum (1 is enough if all chains of that combination are the same) -> MS/MS confirmed.
- Chains: the chain combination (e.g. PC 16:0_18:1) is reported only if fragments for all its chains are present; otherwise species level (PC 34:1).
- Two or more candidates at the same level: all are reported.
- A fragment counts only if it is above the intensity threshold.
- Species of classes or adducts without rules are reported by accurate mass only ("no rules").

## Lyso classes
| Entry | Rule | Source |
|---|---|---|
| LPC [M+H]+ | m/z 184 (diagnostic), m/z 125 and m/z 104 choline (supporting) | Chughtai 2013 (m/z 184, 125); Sugasini 2017 (m/z 104) |
| LPC [M+Na]+, [M+K]+ | m/z 147 / m/z 163 (diagnostic) | Chughtai 2013 |
| LPC [M+Li]+ | NL 59 (diagnostic), NL 183 and NL 189 (supporting) | Hsu & Turk 2003 |
| LPE [M+H]+ | NL 141 (diagnostic) | Ivanova 2007 |
| LPC [M+HCOO]- | NL 60 (diagnostic) | Houjou 2004 (choline lipids as [M+HCOO]- lose HCOO + CH3) |
| LPE [M-H]- | m/z 196 (diagnostic) | Okudaira 2014 |
| LPI [M-H]- | m/z 241 (diagnostic), 259 and 223 (supporting) | Ivanova 2007 |
| LPS [M-H]- | NL 87 (diagnostic) | Ivanova 2007 |
| LPG, LPA | not included | no head-group rule found (Okudaira 2014 uses the FA anion for LPG); accurate mass only |
| Lyso alkali adducts of PE, PS, PI, PG | not included | no source found |

## Not included (not confirmed in the literature)
| Entry | Rule | Reason |
|---|---|---|
| PC [M+K]+ | loss of phosphocholine (183) | not reported for the K adduct |
| PE [M+Na]+ | loss of phosphoethanolamine (141) | not reported for the Na adduct |
| PE [M+K]+ | whole entry | Godzien 2015 did not see a K adduct of PE |
| SM [M+H]+ | loss of H2O; loss of H2O + N(CH3)3 | Hsu & Turk 2000 and Godzien 2015 report only m/z 184 |
| SM [M+K]+ | loss of phosphocholine (183) | not reported for the K adduct |
| PS [M+H]+ | loss of 185 + FA ketene | chain ions weak in MS2 (Hsu & Turk 2005) |
| TG [M+K]+ | whole entry | no source found |
| TG [M+NH4]+ | loss of NH3 | no separate NH3-loss ion (Hsu & Turk 2010) |
| DG [M+NH4]+ | loss of NH3; loss of NH3 + H2O | not reported (Gao 2016 reports only [M-FA]+) |
| DG [M+Li]+ | whole entry | source (Wang 2014) uses derivatized DG |
| MG [M+H]+ | loss of H2O | not reported (Gao 2016 reports only the acylium ion) |
| PC [M+CH3COO]- | m/z 224 | reported only for the formate adduct |
| CL [M-2H]2- | whole entry | FA anion rule not reported; main CL fragments are PA-type ions |
| PGI2 | whole entry | conflict between sources (Dickinson & Murphy 2002); PGI2 is not stable |
| [M-CH3]- of PC | adduct not used | not common in our data |

Classes that are not covered now: CL, PGI2, PE and TG with K. They can be added when a source confirms them.

## Notes on the sources
- On tissue, DG fragments are often not seen at normal collision energy (Lanekoff 2016 SI).
- LPC [M+H]+: the choline ion m/z 104 is present for sn-1 acyl LPC and absent for sn-2 acyl LPC (Sugasini 2017). Thus it is only a supporting fragment.
- Prostaglandin losses: m/z 271 = loss of (2 H2O + CO2), m/z 189 = loss of C7H14O4.

## Per-entry table

Each entry: an example molecule (not in the rule file) and the calculated m/z of each fragment for that molecule. Use these values to check the rules.

## Positive mode

### PC [M+H]+ — example PC 16:0/16:0
Precursor m/z (calc.) **734.5694**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| phosphocholine | ion | diagnostic | 184.0733 |
| cyclophosphane | ion | supporting | 124.9998 |
| loss of phosphocholine | nl | supporting | 551.5034 |
| loss of N(CH3)3 | nl | supporting | 675.4959 |
| loss of FA 16:0 (acid) | nl (chain) | chain | 478.3292 |
| loss of FA 16:0 (ketene) | nl (chain) | chain | 496.3398 |

### PC [M+Li]+ — example PC 16:0/16:0
Precursor m/z (calc.) **740.5776**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of N(CH3)3 | nl | supporting | 681.5041 |
| loss of phosphocholine | nl | diagnostic | 557.5116 |
| loss of lithium cholinephosphate | nl | diagnostic | 551.5034 |
| lithiated cyclophosphane | ion | supporting | 131.0080 |

### PC [M+Na]+ — example PC 16:0/16:0
Precursor m/z (calc.) **756.5514**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of N(CH3)3 | nl | diagnostic | 697.4779 |
| loss of phosphocholine | nl | supporting | 573.4853 |
| loss of sodium cholinephosphate | nl | diagnostic | 551.5034 |
| sodiated cyclophosphane | ion | diagnostic | 146.9818 |

### PC [M+K]+ — example PC 16:0/16:0
Precursor m/z (calc.) **772.5253**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of N(CH3)3 | nl | supporting | 713.4518 |
| loss of potassium cholinephosphate | nl | diagnostic | 551.5034 |
| potassiated cyclophosphane | ion | diagnostic | 162.9557 |

### PE [M+H]+ — example PE 18:0/20:4
Precursor m/z (calc.) **768.5538**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of phosphoethanolamine | nl | diagnostic | 627.5347 |
| acylium FA 18:0 | ion (chain) | chain | 267.2682 |
| acylium FA 20:4 | ion (chain) | chain | 287.2369 |

### PE [M+Li]+ — example PE 16:0/20:4
Precursor m/z (calc.) **746.5307**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of aziridine | nl | supporting | 703.4885 |
| loss of phosphoethanolamine | nl | diagnostic | 605.5116 |
| loss of lithium phosphoethanolamine | nl | diagnostic | 599.5034 |
| [phosphoethanolamine + Li]+ | ion | supporting | 148.0346 |
| loss of FA 16:0 (acid) + C2H5N | nl (chain) | chain | 447.2482 |
| loss of FA 20:4 (acid) + C2H5N | nl (chain) | chain | 399.2482 |

### PE [M+Na]+ — example PE 16:0/18:2
Precursor m/z (calc.) **738.5044**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of aziridine | nl | diagnostic | 695.4622 |
| [phosphoethanolamine + Na]+ | ion | diagnostic | 164.0083 |

### SM [M+H]+ — example SM d18:1/18:0
Precursor m/z (calc.) **731.6062**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| phosphocholine | ion | diagnostic | 184.0733 |

### SM [M+Li]+ — example SM d18:1/16:0
Precursor m/z (calc.) **709.5830**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of N(CH3)3 | nl | supporting | 650.5095 |
| loss of phosphocholine | nl | diagnostic | 526.5170 |
| loss of lithium cholinephosphate | nl | diagnostic | 520.5088 |
| lithiated cyclophosphane | ion | supporting | 131.0080 |

### SM [M+Na]+ — example SM d18:1/16:0
Precursor m/z (calc.) **725.5568**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of N(CH3)3 | nl | diagnostic | 666.4833 |
| loss of phosphocholine | nl | supporting | 542.4908 |
| loss of sodium cholinephosphate | nl | diagnostic | 520.5088 |
| sodiated cyclophosphane | ion | diagnostic | 146.9818 |

### SM [M+K]+ — example SM d18:1/16:0
Precursor m/z (calc.) **741.5307**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of N(CH3)3 | nl | supporting | 682.4572 |
| loss of potassium cholinephosphate | nl | diagnostic | 520.5088 |
| potassiated cyclophosphane | ion | diagnostic | 162.9557 |

### PS [M+H]+ — example PS 16:0/18:1
Precursor m/z (calc.) **762.5280**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of phosphoserine | nl | diagnostic | 577.5190 |

### TG [M+Na]+ — example TG 18:0/18:0/18:0
Precursor m/z (calc.) **913.8195**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of FA 18:0 (acid) | nl (chain) | chain | 629.5479 |
| loss of FA 18:0 (Na salt) | nl (chain) | chain | 607.5660 |

### TG [M+NH4]+ — example TG 18:0/18:0/18:0
Precursor m/z (calc.) **908.8641**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of FA 18:0 (acid) + NH3 | nl (chain) | chain | 607.5660 |

### TG [M+Li]+ — example TG 16:0/18:0/18:1
Precursor m/z (calc.) **867.7987**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of FA 16:0 (acid) | nl (chain) | chain | 611.5585 |
| loss of FA 18:0 (acid) | nl (chain) | chain | 583.5272 |
| loss of FA 18:1 (acid) | nl (chain) | chain | 585.5429 |
| loss of FA 16:0 (Li salt) | nl (chain) | chain | 605.5503 |
| loss of FA 18:0 (Li salt) | nl (chain) | chain | 577.5190 |
| loss of FA 18:1 (Li salt) | nl (chain) | chain | 579.5347 |

### DG [M+NH4]+ — example DG 16:0/18:1
Precursor m/z (calc.) **612.5562**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of FA 16:0 (acid) + NH3 | nl (chain) | chain | 339.2894 |
| loss of FA 18:1 (acid) + NH3 | nl (chain) | chain | 313.2737 |

### DG [M+Na]+ — example DG 14:0/14:0
Precursor m/z (calc.) **535.4333**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| [FA 14:0+Na]+ | ion (chain) | chain | 251.1982 |

### DG [M+K]+ — example DG 14:0/14:0
Precursor m/z (calc.) **551.4072**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| [FA 14:0+K]+ | ion (chain) | chain | 267.1721 |

### MG [M+H]+ — example MG 17:1
Precursor m/z (calc.) **343.2843**. Class rule: chains.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| acylium FA 17:1 | ion (chain) | chain | 251.2369 |

## Negative mode

### PC [M+HCOO]- — example PC 16:0/18:2
Precursor m/z (calc.) **802.5604**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of HCOO + CH3 ([M-CH3]-) | nl | diagnostic | 742.5392 |
| demethylated phosphocholine | ion | supporting | 168.0431 |
| demethylated glycerophosphocholine - H2O | ion | supporting | 224.0693 |
| FA 16:0 anion | ion (chain) | chain | 255.2330 |
| FA 18:2 anion | ion (chain) | chain | 279.2330 |

### PC [M+CH3COO]- — example PC 17:0/20:4
Precursor m/z (calc.) **854.5917**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of CH3COO + CH3 ([M-CH3]-) | nl | diagnostic | 780.5549 |
| demethylated phosphocholine | ion | supporting | 168.0431 |
| FA 17:0 anion | ion (chain) | chain | 269.2486 |
| FA 20:4 anion | ion (chain) | chain | 303.2330 |

### PE [M-H]- — example PE 16:0/18:1
Precursor m/z (calc.) **716.5236**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| phosphoethanolamine - H | ion | diagnostic | 140.0118 |
| glycerophosphoethanolamine - H2O - H | ion | diagnostic | 196.0380 |
| FA 16:0 anion | ion (chain) | chain | 255.2330 |
| FA 18:1 anion | ion (chain) | chain | 281.2486 |
| loss of FA 16:0 (ketene) | nl (chain) | chain | 478.2939 |
| loss of FA 18:1 (ketene) | nl (chain) | chain | 452.2783 |

### PI [M-H]- — example PI 18:0/20:4
Precursor m/z (calc.) **885.5499**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| inositol phosphate - H2O | ion | diagnostic | 241.0119 |
| inositol phosphate | ion | supporting | 259.0224 |
| inositol phosphate - 2H2O | ion | supporting | 223.0013 |
| FA 18:0 anion | ion (chain) | chain | 283.2643 |
| FA 20:4 anion | ion (chain) | chain | 303.2330 |
| loss of FA 18:0 (acid) | nl (chain) | chain | 601.2783 |
| loss of FA 20:4 (acid) | nl (chain) | chain | 581.3096 |

### PS [M-H]- — example PS 16:0/18:1
Precursor m/z (calc.) **760.5134**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of serine | nl | diagnostic | 673.4814 |
| FA 16:0 anion | ion (chain) | chain | 255.2330 |
| FA 18:1 anion | ion (chain) | chain | 281.2486 |
| loss of FA 16:0 (acid) + C3H5NO2 | nl (chain) | chain | 417.2411 |
| loss of FA 18:1 (acid) + C3H5NO2 | nl (chain) | chain | 391.2255 |

### PG [M-H]- — example PG 16:0/18:1
Precursor m/z (calc.) **747.5182**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| glycerophosphate | ion | diagnostic | 171.0064 |
| glycerophosphate - H2O | ion | supporting | 152.9958 |
| FA 16:0 anion | ion (chain) | chain | 255.2330 |
| FA 18:1 anion | ion (chain) | chain | 281.2486 |
| loss of FA 16:0 (acid) | nl (chain) | chain | 491.2779 |
| loss of FA 18:1 (acid) | nl (chain) | chain | 465.2623 |

### PGE2 / PGD2 / PGH2 [M-H]- — example PGE2 / PGD2 / PGH2
Precursor m/z (calc.) **351.2177**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of H2O | nl | supporting | 333.2071 |
| loss of 2H2O + CO2 | nl | supporting | 271.2067 |
| loss of C7H14O4 | nl | diagnostic | 189.1285 |

### PGF2a [M-H]- — example PGF2a
Precursor m/z (calc.) **353.2333**. Class rule: head group.

| fragment | type | role | calc m/z (example) |
|---|---|---|---|
| loss of C2H4O | nl | supporting | 309.2071 |
| C12H17O2- | ion | diagnostic | 193.1234 |
