---
title: Sample size estimate in T-test
author: Junli Zhang
date: '2017-12-16'
slug: sample-size-estimate-in-t-test
categories:
  - Research
tags:
  - math
  - lessons
---

I planted 16 F2 plants to test whether a marker is significant in this population. I get about 4 plants of AA and 4 plants of BB, and I did the T-test between AA and BB group, but I found the pvalue is only about 0.08. Should I reject the null hypothesis that AA and BB are the same? I probably will, but this experiment has two few samples, so I do not have enough statistical power to reject the null hypothesis.

After a search, I found we could estimate the sample size we need in order to get a significance level, say 0.001. We just need to provide the desired significance level, the sample mean difference between the two groups (based on our first experiment or some other prior), the standard deviation of the population, and the desired power (`$1 - \beta$`, the possibility it is real, `$\beta$` is type II error). [Here is some detailed information](https://stats.idre.ucla.edu/other/gpower/power-analysis-for-two-group-independent-sample-t-test/).

![t-test](/images/20171217_ttest.jpg)

There is also [a good website for sample estimate](https://www.stat.ubc.ca/~rollin/stats/ssize/n2.html). Here I just remade the wheel for my practice.

<script src="/libs/pvalue_to_sample_size_ttest.js"></script>
<div id="sample_size">

<p>Enter parameters below to estimate required sample size:</p>

<FORM>
<pre>
Mean Difference                  : <INPUT TYPE="text" NAME="MD" Value="2.0" SIZE=15>
Alpha (default 0.05, right tail) : <INPUT TYPE="text" NAME="alpha" Value="0.05" SIZE=15>
Standard Deviation               : <INPUT TYPE="text" NAME="stdev" Value="1.0" SIZE=15>
Desired power (default 0.80)     : <INPUT TYPE="text" NAME="power" Value="0.8" SIZE=15>
</pre>
The sample size (for each sample separately): <INPUT TYPE="text" NAME="result" SIZE=20> <INPUT TYPE="button" VALUE="Calculate" ONCLICK="compute(this.form)">
</FORM>
</div>

For example, in the experiment I mentioned, the group mean difference between AA and BB is about 2, standard deviation is about 1.2, and I want significant level `$\alpha = 0.001$` (right tail, because I want AA - BB = 2), so the required sample size for each group is about 11, therefor I need to plant about `11 x 4 = 44` F2 plants.
