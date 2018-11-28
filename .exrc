augroup exrc
    autocmd!
    autocmd Bufread examples/* let b:start = 'python -m '.substitute(expand('%:r'), '\/', '.', 'g')
    autocmd Bufread doc/* set makeprg=make\ documentation
augroup END

