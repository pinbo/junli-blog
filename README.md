# junli-blog

My blog to share experiences, knowledge and ideas.

It is hosted here:

https://junli.netlify.app/

# build website with hugo
```
# since v0.60.0, hugo render was switched from Blackfriday to goldmark
# options(blogdown.hugo.version = "0.59.1") 
library(blogdown)
install_hugo(version = "0.59.1")
serve_site() # hugo server --bind 127.0.0.1 -p 4321 -D -F
stop_server()
```
