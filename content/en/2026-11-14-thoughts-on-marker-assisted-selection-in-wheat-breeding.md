---
title: Thoughts on Marker-assisted Selection in Wheat Breeding
author: Junli Zhang
date: '2026-11-14'
slug: thoughts-on-marker-assisted-selection-in-breeding
categories:
  - thoughts
tags:
  - Marker-assisted selection
  - breeding
mermaid: true
---

Marker-Assisted Selection (MAS) is a great tool to select lines with specific genes. I have been using it in my genetic studies: 
1) to select recombinant lines between two genetic markers; 
2) to select BCnF1 lines for continuing backcrossing. 

But I am not clear on how to use it in a real breeding program. I thought I knew, but I found I did not after working alongside wheat breeders in my current position. Initially, I considered using MAS to select F2 lines with specific genes. However, in practice, this approach faces significant challenges:
1. Breeders typically produce hundreds of F1 crosses annually. Even if only 100 F2 plants are produced per cross, the volume of DNA extraction and genotyping required would be immense.
2. The probability (`$1/2^n$`, `$n$` is the number of unlinked markers) is very low when selecting an F2 line that carries multiple desired markers simultaneously. Consequently, breeders would have to discard many lines, inadvertently eliminating genetic variation for selection.

Although wheat curl mite and greenbug are a big problem in the Texas High Plains, not all released varieties have the *Cmc4* and *Gb3* genes. I guess this is because wheat curl mite and greenbug do not occur every year during the selection process, so there is no selection pressure. Or maybe there are just no specific combinations for breeders to select for.

This August, while attending the Texas Small Grains Workers’ meeting, I had the opportunity to see how breeders decide which lines to proceed with and which to discard. Some lines were discarded because they did not perform well in areas with high pest pressure (such as Hessian fly in South and East Texas this year), although they performed well when pest pressure was low. These lines would perform well if we can add the required genes. ***I think this is the stage we could fully use MAS, more specifically marker-assisted backcrossing (MAB)***:
Before their final release, we could use MAB to fix the problematic lines on disease resistance or end-use quality. With speed breeding, this adds only 1 to 2 more years’ work to the breeding process, but saves 5 to 6 years by avoiding a new breeding cycle. In this context, our goal is to simply introduce necessary disease-resistance genes (e.g., Hessian fly resistance) into an otherwise superior candidate. Since these candidates already possess several years of preliminary trial data, they are well on their way to field trials after incorporating the new genes.

To maximize efficiency, the ideal approach is to prepare a single crossing donor that carries the comprehensive set of required genes for the target region. For example, for the Texas High Plains, the donor should possess resistance to wheat curl mite (*Cmc4*), wheat streak mosaic virus (*wsm2* and *wsm3*), leaf rust (*Lr21* and *Lr46*), stripe rust (*Yr5* and *Yr15*), Hessian fly (*H13*, *H21*, *H33*), and greenbug (*Gb3*). Ideally, the donor should carry two genes per resistance trait to ensure resistance durability. Once we create a donor with all 10 genes mentioned above, we can perform a single cross and select for all traits simultaneously in the offspring. However, as the number of genes we wish to stack increases, the number of pollinations required to generate enough seed for the next selection round becomes large. If we target 4 genes, it is still feasible (1/16 of the offspring will have all 4 genes). As this practice continues, most of the crossing pool should have those genes, so it will hopefully be easier in the future.

**Wheat breeding process in Texas High Plains**

{{< mermaid >}}
flowchart TB
    %% breeding years
    A(["Y1: Crossing block"])
    B(["Y2: F1 (200)"])
    C(["Y3-5: F2-F4 bulk plots (200-450)"])
    D(["Y6: F4:5 Head-Rows (20,000)"])
    E(["Y7: F4:6 Observation Yield Trials (600)"])
    F(["Y8: Preliminary Yield Trials (200)"])
    G(["Y9: Advanced Yield Trials (120)"])
    H(["Y10: Texas Elite Yield Trial (TXE, 15)"])
    I(["Y11: SRPN, TXE, and increase (6)"])
    J(["Y12: SRPN, TXE, UVT, increase, and WQC (2)"])
    K(["Y13: Release (1)"]):::decision
    L("Marker-assisted backcrossing to introduce missing genes")
    %% Phase 1: Crossing & Early Generation
    subgraph P1 [Phase 1: Crossing & Segregation]
        direction LR
        A --> B --> C --> D
    end

    %% Yield Trials
    subgraph P2 [Phase 2]
        direction TB
        subgraph YT [Yield trials]
            E --> F --> G --> H
        end
        L --> YT --> L
    end

    subgraph P3 [Phase 3: more test and release]
        direction LR
        I --> J --> K
    end
    P1 --> P2 --> P3
{{< /mermaid >}}
