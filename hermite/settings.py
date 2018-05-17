settings = {
        'cache': False,
        'cachedir': '/tmp',
        'tensorize': True,
        'trails': False,
        'debug': False
        }


def get(name, **kwargs):
    if name in kwargs:
        return kwargs['tensorize']
    return settings[name]
