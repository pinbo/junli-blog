---
title: Two good tutorials for creating Hugo websites through GitHub pages and Netlify
author: Junli Zhang
date: '2017-11-18'
slug: two-good-tutorials-for-creating-hugo-websites-through-github-pages-and-netlify
categories:
  - Hugo
tags:
  - Tips
---

I found two great tutorials using Hugo to build your own website and then publish on either GitHub or Netlify.

1. Publish on github with [Hugo](https://gohugo.io/) directly: https://youtu.be/3wkR8GyDODs

1. Publish on [Netlify](https://www.netlify.com/) with [Blogdown](https://bookdown.org/yihui/blogdown/):   https://alison.rbind.io/post/up-and-running-with-blogdown/

One update for the 2nd tutorial is that for build command:
https://www.netlify.com/docs/continuous-deployment/#common-configuration-directives

```
For Hugo hosting, the build command hugo will build and deploy with the version 0.17 of Hugo. For versions 0.13, 0.14, 0.15, 0.16, 0.17, 0.18 and 0.19, you can specify a specific Hugo release like this: hugo_0.15. For version 0.20 and above, use the regular hugo command and create a Build Environment Variable called HUGO_VERSION and set it to the version of your choice. (More details can be found on our blog.)
```

For me, I am using hugo v0.29, so I cannot use `hugo_0.29` as the command'; instead I need to just use `hugo` as the command, and create a **Build Environment Variable** called `HUGO_VERSION` and set the value `0.29`. This can be done in the setting page.
