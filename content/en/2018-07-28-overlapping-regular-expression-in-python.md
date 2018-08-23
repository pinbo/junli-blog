---
title: Overlapping regular expression in python
author: Junli Zhang
date: '2018-07-28'
slug: overlapping-regular-expression-in-python
categories: [programming]
tags:
  - python
---

I did not realize that by default python only find patterns in a non-overlapping way. For example, we want to find "AA" in "AAAA", it will give you positions 0 and 2, eg the first two AA and the last two AA. However, the middle two AA is also a match. To get all the overlapping matches, we need "lookahead assertion". From the python regex doc:

> (?=...)  
> Matches if ... matches next, but doesn’t consume any of the string. This is called a lookahead assertion. For example, Isaac (?=Asimov) will match 'Isaac ' only if it’s followed by 'Asimov'.

A great example is in this post:  https://stackoverflow.com/questions/5616822/python-regex-find-all-overlapping-matches

So to find all the starting positions of the pattern, we need to put the pattern to `(?=...)`:

``` python
import re
ss = "AAAA"
for m in re.finditer('(?=(AA))', ss):
    print m.start(), m.end(), m.group(1)
```

`m.start()` and `m.end()` are the same, because actually our pattern is an empty string "", `m.group(1)` will give us the string matched the pattern "AA".

`re.findall(r'(?=(AA))', ss)` and `re.findall(r'(?=AA)', ss)` give different results: the first one gives the matched "AA", but the 2nd one gives the matched "", but the m.start() is the same. I have not figured out why.
