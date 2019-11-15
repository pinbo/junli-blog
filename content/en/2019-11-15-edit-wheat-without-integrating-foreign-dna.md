---
title: Edit wheat without integrating foreign DNA
author: Junli Zhang
date: '2019-11-15'
slug: edit-wheat-without-integrating-foreign-dna
categories:
  - Ideas
tags:
  - transgenic
  - wheat
---

Now to edit wheat (possibly most other plants), we need to use Agrobacteria or Gene gan to integrate vectors with CAS9 and selection markers into wheat's genome, then select plants with editing events but without the vectors. Usually, the insertion of vectors happens in multiple places of the genome (multi-copy), so we need a big number of F2 to select the edited but non-GMO plants. We still use markers to see whether there are vector sequences in the genome, so we cannot be 100% sure that not a little foreign DNA present in the wheat genome unless we do the whole genome sequencing. All these things let me think about how to edit the plant genome without insertion of the vector sequences (transient method). There are already some transient methods to edit the wheat genome. Here are what I know so far.

1. Protoplast transformation with small vectors. This is actually the way how we test the editing efficiency of the guide RNAs. We can see editing events easily in protoplasts, but there is no way to regenerate a plant from the protoplast for wheat.

2. Pollinate wheat with transgenic corn. The paper ["One-step genome editing of elite crop germplasm during haploid induction"](https://www.nature.com/articles/s41587-019-0038-x) proved this is possible. We use this method to make double haploid wheat because corn pollen can enter the wheat embryo but cannot merge with the mother cell. One paper has proved that editing events can happen, but the seeds are haploid. This method is good but we need to transform corn first and the editing efficiency is still low.

3. Introduce CAS9 + sgRNA construct into immature wheat embryos by particle bombardment. Particle bombardment, however, may break the genome and still introduce the foreign genes into the wheat genome, as shown by [Zhang et al. 2016](https://www.nature.com/articles/ncomms12617) that the CRISPR/Cas9 DNA construct was found to be present about 50% of the T0 plants.

4. Deliver Cas9 and sgRNA *in vitro* transcripts (IVTs) into the immature embryos. This method will eliminate the possibility of integrating vector DNA into the genome, but it has very low efficiency (1.1%) ([Zhang et al. 2016](https://www.nature.com/articles/ncomms12617)), so thousands of plants need to be screened to get desirable mutations. Besides, you also need to prepare IVTs.

5. Directly deliver [CRISPR/Cas9 ribonucleoprotein (RNP) complexes](https://www.nature.com/articles/ncomms14261) into the immature embryo cells of the bread wheat by particle bombardment. It seems it can achieve 4% to 5% of edited embryos. It looks like a promising method, but not widely used yet.

6. Use nanotubes to deliver vectors to wheat. The paper ["High aspect ratio nanomaterials enable delivery of functional genetic material without DNA integration in mature plants"](https://www.nature.com/articles/s41565-019-0382-5) reported successful cases in wheat, but it is not repeated by others yet.

7. The method I am thinking about. Edit the Agrobacterial T-vector so the vector can be sent to the nucleus but cannot integrate into the genome. I think this is possible, but I do not have time to investigate.