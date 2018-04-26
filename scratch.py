def natural_bissect(func, x1=0, x2=1000):
    f1, f2 = func(x1), func(x2)
    if f1 is 0:
        return x1
    elif f2 is 0:
        return x2
    assert f1*f2 < 0
    x3 = (x1+x2)//2
    f3 = func(x3)
    if f1*f3 <= 0:
        return natural_bissect(func, x1=x1, x2=x3)
    else:
        return natural_bissect(func, x1=x3, x2=x2)

print(natural_bissect(lambda x: (x-2)*(x-10)*(x-25)))
