---
title: Try using remark.js markdown style
author: Junli Zhang
date: '2018-03-25'
slug: try-using-remark-js-markdown-style
categories:
  - Linux
tags:
  - blogdown
  - remark
---

Remark.js for presentation has many easy-to-use features, especially its powerful macros. So I want to use remark.js to convert my markdown to HTML in blogdown (or Hugo) to avoid remembering different formatting styles. So I made a new template to directly use `remark.js` to create the final HTML page.

`![:color green](This should be green color!)` will show ![:color green](This should be green color!)

OR

`I have .red[red color]` will show "I have .red[red color]".

It cannot use blocks right now.

### Math is supported too
When \\(a \ne 0\\), there are two solutions to \\(ax^2 + bx + c = 0\\) and they are
  $x = {-b \pm \sqrt{b^2-4ac} \over 2a}.$


