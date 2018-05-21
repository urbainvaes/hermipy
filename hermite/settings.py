settings = {
        'cache': False,
        'cachedir': '.cache',
        'tensorize': True,
        'trails': False,
        'debug': False
        }


def get(name, **kwargs):
    if name in kwargs:
        return kwargs[name]
    return settings[name]
