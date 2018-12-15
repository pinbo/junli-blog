---
title: Nvidia driver for Ubuntu 18.10 on Thinkpad T480
author: Junli Zhang
date: '2018-12-14'
slug: nvidia-driver-for-ubuntu-18-10-on-thinkpad-t480
categories:
  - Linux
tags:
  - nvidia
  - T480
---

I installed Ubuntu 18.10 on my new Thinkpad T480. However, the system froze every time after I logged in. It turned out to be the problem of Nvidia driver. After its freezing, I also could not use **Ctrl+Alt+F1** to use the shell to reboot. But later I realized that I need to turn on `FnLock` to use the F1 to F12 keys. Finally, [this post](https://askubuntu.com/questions/1044813/ubuntu-18-04-and-nvidia-stuck-after-boot) solved my problem. The steps are below:

1. Turn on the laptop and press `Shift` or `Esc` during boot to show the GRUB menu
2. Press `e` and edit the line starting with `linux` adding `nouveau.modeset=0` at the end of the line
3. Press `F10` to boot
4. Once logged in successfully, go to **Activities** (top left corner) **Drivers** and install suggested Nvidia driver.

Now it should work without any problems. I did not check whether the GRUB menu still has my edits.