augroup exrc
    autocmd!
    autocmd Bufread examples/* let b:start = 'python -m '.substitute(expand('%:r'), '\/', '.', 'g')
augroup END
