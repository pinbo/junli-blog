---
title: Make slides work in blogdown
author: Junli Zhang
date: '2017-11-16'
slug: make-slides-work-in-blogdown
categories:
  - Hugo
tags:
  - Tips
---

I made some slides introducing primer design tips. I want to put it in my blog; however, blogdown cannot render [xaringan](https://github.com/yihui/xaringan) Rmd files if I put them in the post folder. I searched a little bit and found one solution: put the raw md file content in your remark.js slides html file. Here is how I do it. I am using the hugo-lithium-theme as an example.

- I made the slides using [xaringan](https://github.com/yihui/xaringan);

- Based on the slide html file, I made a [**head_remarkjs.html**](/examples/head_remarkjs.html) and a  [**foot_remarkjs.html**](/examples/foot_remarkjs.html) in folder `themes/hugo-lithium-theme/layouts/partials`

- Modified the single.html file located in `themes/hugo-lithium-theme/layouts/_default` by adding an ifelse condition to see if the markdown file has tag "slides". If yes, it will make a single file using `head_remarkjs.html + raw content of the md file + foot_remarkjs.html`; else it will use the default single.html layout.

Here is the header I added to the top of the single.html file:
```html
{{ if in .Params.tags "slides" }}

{{/* Slides mode */}}
{{ partial "head_remarkjs.html" . }}
	<body>
	<textarea id="source">
{{ .RawContent }}
	</textarea>
{{ partial "foot_remarkjs.html" . }}

{{else}}

```

Make sure your slide md file has tag "slides". You can put any other images in a folder named "images" (or other names if you like) in the **static/** folder, and use `[my image](/image/myimage.png)` to refer it in the md file.

Later I found Yihui already made an introduction how to using [xaringan](https://github.com/yihui/xaringan) in [blogdown](https://github.com/rstudio/blogdown). Here is the website:
https://github.com/yihui/blogdown-static

This method is much easier, and you just need to:

1. Create a **R/** folder and add an R script `build.R` with one line of code `blogdown::build_dir('static')` in the  **R/** folder

1. Put your Rmd files under the **static/** folder.

That is it! The one thing I am not satisfied with this method is that the slides will not be shown on my blog list. So I made another file `slides.md` in the `content/` folder along with the `about.md` to list out all my slides:

```md
---
title: "Slides"
date: "2017-11-16"
---

Here a list of slides I made.

- [slide 1](/slides/slide1.html)

- [slide 2](/slides/slide2.html)
```

And add 1 menu for the slides on my homepage by adding:
```
[[menu.main]]
    name = "Slides"
    url = "/slides/"
```
in the file `config.toml` with other `[[menu.main]]`.

I am just a beginner to use blogdown and hugo. I am sure there will be better ways. Please do not hesitate to share yours.
