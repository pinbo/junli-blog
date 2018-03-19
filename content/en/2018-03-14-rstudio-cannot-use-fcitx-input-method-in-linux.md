---
title: Rstudio cannot use fcitx input method in Linux
author: Junli Zhang
date: '2018-03-14'
slug: rstudio-cannot-use-fcitx-input-method-in-linux
categories:
  - Linux
tags:
  - setting
---

I am currently using Rstudio Version 1.1.419, which cannot use Chinese input method with fcitx. The [official solution](https://support.rstudio.com/hc/en-us/articles/205605748-Using-RStudio-0-99-with-Fctix-on-Linux) does not work anymore. Finally I found this [blog](https://bbs.deepin.org/forum.php?mod=viewthread&tid=149730) that fixed the problem. Here is my summary on Debian stretch:

```{bash}
wget http://ikuya.info/tmp/fcitx-qt5-rstudio.tar.gz 
tar xf fcitx-qt5-rstudio.tar.gz 
cd fcitx-qt5-rstudio 
sudo dpkg -i libfcitx-qt5-1-rstudio_1.0.5-1_amd64.deb

# installation of fcitx-frontend-qt5-rstudio_1.0.5-1_amd64.deb has a problem
# here we just get the library
dpkg -X fcitx-frontend-qt5-rstudio_1.0.5-1_amd64.deb extract
sudo cp extract/usr/lib/rstudio/bin/plugins/platforminputcontexts/libfcitxplatforminputcontextplugin.so /usr/lib/rstudio/bin/plugins/platforminputcontexts/
```

Hope this helps and hope Rstudio can solve the problem in the next version!