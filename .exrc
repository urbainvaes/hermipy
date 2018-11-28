augroup exrc
    autocmd!
    autocmd Bufread examples/* let b:start = 'python -m '.substitute(expand('%:r'), '\/', '.', 'g')
    autocmd Bufread docs/* set makeprg=make\ doc
augroup END

