---
title: Setup fcitx5 for Chinese input in WPS and WeChat
author: Junli Zhang
date: '2025-01-01'
slug: setup-fcitx5-for-wps-and-wechat
categories:
  - linux
  - method
tags:
  - kde
  - fcitx
  - wps
  - wechat
  - chinese
---

Recently, I installed [Linux CachyOS](https://cachyos.org/) v6.12.7 with KDE Plasma v6.2.4. I write a lot of Chinese, so I installed [fcitx5](https://fcitx-im.org/wiki/Fcitx_5) with Chinese input methods:

``` sh
sudo pacman -S fcitx5-im fcitx5-chinese-addons
```

The Chinese input methods work in Google Chrome, Firefox, and other applications, but not for WPS Office and WeChat. I checked the [wiki page of fcitx5](https://wiki.archlinux.org/title/Fcitx5), and used the command `fcitx5-diagnose` to see whether I can get some clues. I made a '.xprofile' file and put the commands below:

``` sh
export XMODIFIERS=@im=fcitx
export QT_IM_MODULE=fcitx
export GTK_IM_MODULE=fcitx
```

However, I still could not activate Chinese input methods in WPS and WeChat. I googled and found that I need to add some "environment variables" when running these two applications. I tried in terminal and it works. Later I found an easier method to add "environment variables" for any programs in the START menu (right click the entry and choose `Edit Application`).
![modify menu entry in KDE start menu](/images/20250101_KDE_start_menu.png)

Just add `XMODIFIERS=@im=fcitx` in the "Environment Variables" field for WeChat and `QT_IM_MODULE=fcitx` for WPS (if writer, spreadsheet, and presentation are separate apps, modify all of them).

![modiy wechat environment variable](/images/20250101_wechat_desktop_change.png)

![modiy wps environment variable](/images/20250101_wps_desktop_change.png)