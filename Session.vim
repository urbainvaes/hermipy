let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/scratch/integrator
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +1 libhermite/src/hermite/transform.cpp
badd +107 mckean-vlasov.py
badd +1 spectral.py
badd +1749 term://.//27921:ipython
badd +1 /tmp/home-urbain-scratch-integrator-Fokker-Planck.ipynb
badd +0 ~/personal/vim/.vim/mySnippets/ipynb.snippets
argglobal
silent! argdel *
edit libhermite/src/hermite/transform.cpp
set splitbelow splitright
set nosplitbelow
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=marker
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 37 - ((27 * winheight(0) + 21) / 43)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
37
normal! 0
lcd ~/scratch/integrator
tabedit ~/scratch/integrator/mckean-vlasov.py
set splitbelow splitright
set nosplitbelow
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=marker
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
1
normal! zo
14
normal! zo
146
normal! zo
let s:l = 8 - ((7 * winheight(0) + 21) / 42)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
8
normal! 018|
lcd ~/scratch/integrator
tabedit ~/scratch/integrator/spectral.py
set splitbelow splitright
set nosplitbelow
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=marker
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 21) / 42)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
lcd ~/scratch/integrator
tabedit ~/scratch/integrator/mckean-vlasov.py
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
set nosplitbelow
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 63 + 63) / 127)
exe '2resize ' . ((&lines * 21 + 22) / 45)
exe 'vert 2resize ' . ((&columns * 63 + 63) / 127)
exe '3resize ' . ((&lines * 20 + 22) / 45)
exe 'vert 3resize ' . ((&columns * 63 + 63) / 127)
argglobal
setlocal fdm=marker
setlocal fde=0
setlocal fmr=```{.,```
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 2 - ((1 * winheight(0) + 21) / 42)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 0
lcd ~/scratch/integrator
wincmd w
argglobal
if bufexists('~/personal/vim/.vim/mySnippets/ipynb.snippets') | buffer ~/personal/vim/.vim/mySnippets/ipynb.snippets | else | edit ~/personal/vim/.vim/mySnippets/ipynb.snippets | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr=```{.,```
setlocal fdi=#
setlocal fdl=99
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 2 - ((1 * winheight(0) + 10) / 21)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 010|
lcd ~/scratch/integrator
wincmd w
argglobal
if bufexists('/tmp/home-urbain-scratch-integrator-Fokker-Planck.ipynb') | buffer /tmp/home-urbain-scratch-integrator-Fokker-Planck.ipynb | else | edit /tmp/home-urbain-scratch-integrator-Fokker-Planck.ipynb | endif
setlocal fdm=marker
setlocal fde=0
setlocal fmr=```{.,```
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
1
normal! zo
11
normal! zo
let s:l = 1 - ((0 * winheight(0) + 10) / 20)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
lcd ~/scratch/integrator
wincmd w
3wincmd w
exe 'vert 1resize ' . ((&columns * 63 + 63) / 127)
exe '2resize ' . ((&lines * 21 + 22) / 45)
exe 'vert 2resize ' . ((&columns * 63 + 63) / 127)
exe '3resize ' . ((&lines * 20 + 22) / 45)
exe 'vert 3resize ' . ((&columns * 63 + 63) / 127)
tabnext 4
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
let g:this_session = v:this_session
let g:this_obsession = v:this_session
let g:this_obsession_status = 2
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
