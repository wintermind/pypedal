#coding=utf-8
# dump python dict to ini format file
# Author: limodou (limodou@gmail.com)
# Copyleft GPL
# $Revision: 245 $
# you can see http://wiki.woodpecker.org.cn/moin/Dict4Ini for more details
#
# Updates:
#   2005/12/12
#     Fixed "\" bug in option's value
#   2005/12/09
#     Adding dict() method, then you can change the DictIni object to dict type
#   2005/10/16
#     Saving the order of the items
#     Adding float format
#

__version__ = '0.4'

import sys
import locale
import os.path
import re

r_float = re.compile('\d*\.\d+')
section_delimeter = '/'

class DictNode(object):
    def __init__(self, values, encoding=None, root=None, section=[], orders=[]):
        self._items = values
        self._orders = orders
        self._encoding = encoding
        self._root = root
        self._section = section

    def __getitem__(self, name):
        if self._items.has_key(name):
            value = self._items[name]
            if isinstance(value, dict):
                return DictNode(value, self._encoding, self._root, self._section + [name])
            else:
                return value
        else:
            self._items[name] = {}
            self._root.setorder(self.get_full_keyname(name))
            return DictNode(self._items[name], self._encoding, self._root, self._section + [name])

    def __setitem__(self, name, value):
        if section_delimeter in name:
            sec = name.split(section_delimeter)
            obj = self._items

            _s = self._section[:]
            for i in sec[:-1]:
                _s.append(i)
                if obj.has_key(i):
                    if isinstance(obj[i], dict):
                        obj = obj[i]
                    else:
                        obj[i] = {} #may lost some data
                        obj = obj[i]
                else:
                    obj[i] = {}
                    self._root.setorder(section_delimeter.join(_s))
                    obj = obj[i]
            obj[sec[-1]] = value
            self._root.setorder(section_delimeter.join(_s + [sec[-1]]))
        else:
            self._items[name] = value
            self._root.setorder(self.get_full_keyname(name))

    def __delitem__(self, name):
        if self._items.has_key(name):
            del self._items[name]

    def __repr__(self):
        return repr(self._items)

    def __getattr__(self, name):
        return self.__getitem__(name)

    def __setattr__(self, name, value):
        if name.startswith('_'):
            if name == '_comment':
                self._root._comments[section_delimeter.join(self._section)] = value
            else:
                self.__dict__[name] = value
        else:
            self.__setitem__(name, value)

    def comment(self, name, comment):
        if name:
            self._root._comments[section_delimeter.join(self._section + [name])] = comment
        else:
            self._root._comments[section_delimeter.join(self._section)] = comment

    def __delattr__(self, name):
        if self._items.has_key(name):
            del self._items[name]

    def __str__(self):
        return repr(self._items)

    def __len__(self):
        return len(self._items)

    def has_key(self, name):
        return self._items.has_key(name)

    def items(self):
        return self._items.items()

    def setdefault(self, name, value):
        return self._items.setdefault(name, value)

    def get(self, name, default=None):
        return self._items.get(name, default)

    def keys(self):
        return self._items.keys()

    def values(self):
        return self._items.values()

    def get_full_keyname(self, key):
        return section_delimeter.join(self._section + [key])

class DictIni(DictNode):
    def __init__(self, inifile=None, values=None, encoding=None, commentdelimeter='#'):
        self._items = {}
        self._inifile = inifile
        self._root = self
        self._section = []
        self._commentdelimeter = commentdelimeter
        self._comments = {}
        self._orders = {}
        self._ID = 1
        self._encoding = getdefaultencoding(encoding)
        if values is not None:
            self._items = values

        if self._inifile and os.path.exists(self._inifile):
            self.read(self._inifile, self._encoding)
            
    def setfilename(self, filename):
        self._inifile = filename

    def getfilename(self):
        return self._inifile

    def save(self, inifile=None, encoding=None):
        if inifile is None:
            inifile = self._inifile

        if isinstance(inifile, (str, unicode)):
            f = file(inifile, 'w')
        elif isinstance(inifile, file):
            f = inifile
        else:
            f = inifile

        if not f:
            f = sys.stdout

        if encoding is None:
            encoding = self._encoding

        f.write(self._savedict([], self._items, encoding))
        if isinstance(inifile, (str, unicode)):
            f.close()

    def _savedict(self, section, values, encoding):
        if values:
            buf = []
            default = []
            for key, value in self._getorderitems(values.items()):
                if isinstance(value, dict):
                    sec = section[:]
                    sec.append(key)
                    buf.append(self._savedict(sec, value, encoding))
                else:
                    c = self._comments.get(section_delimeter.join(section + [key]), '')
                    if c:
                        lines = c.splitlines()
                        default.append('\n'.join(['%s %s' % (self._commentdelimeter, x) for x in lines]))

                    default.append("%s = %s" % (key, uni_prt(value, encoding)))
            if default:
                buf.insert(0, '\n'.join(default))
                buf.insert(0, '[%s]' % section_delimeter.join(section))
                c = self._comments.get(section_delimeter.join(section), '')
                if c:
                    lines = c.splitlines()
                    buf.insert(0, '\n'.join(['%s %s' % (self._commentdelimeter, x) for x in lines]))
            return '\n'.join(buf + [''])
        else:
            return ''
        
    def dict(self):
        return self._dict(self)
                
    def _dict(self, v):
        if isinstance(v, tuple):
            return tuple([self._dict(x) for x in v])
        elif isinstance(v, list):
            return [self._dict(x) for x in v]
        elif isinstance(v, (dict, DictNode)):
            d = {}
            for key, value in v.items():
                d[key] = self._dict(value)
            return d
        else:
            return v
        
    def read(self, inifile=None, encoding=None):
        if inifile is None:
            inifile = self._inifile

        if isinstance(inifile, (str, unicode)):
            try:
                f = file(inifile, 'r')
            except:
                return  #may raise Exception is better
        elif isinstance(inifile, file):
            f = inifile
        else:
            f = inifile

        if not f:
            f = sys.stdin

        if encoding is None:
            encoding = self._encoding

        comments = []
        section = ''
        for line in f.readlines():
            line =  line.strip()
            if not line: continue
            if line.startswith(self._commentdelimeter):
                comments.append(line[1:].lstrip())
                continue
            if line.startswith('['):    #section
                section = line[1:-1]
                #if comment then set it
                if comments:
                    self.comment(section, '\n'.join(comments))
                    comments = []
                continue
            key, value = line.split('=', 1)
            key = key.strip()
            value = process_value(value.strip(), encoding)
            if section:
                self.__setitem__(section + section_delimeter + key, value)
                #if comment then set it
                if comments:
                    self.__getitem__(section).comment(key, '\n'.join(comments))
                    comments = []
            else:
                self.__setitem__(key, value)
                #if comment then set it
                if comments:
                    self.comment(key, '\n'.join(comments))
                    comments = []
        if isinstance(inifile, (str, unicode)):
            f.close()

    def setorder(self, key):
        if not self._orders.has_key(key):
            self._orders[key] = self._ID
            self._ID += 1

    def _getorderitems(self, values):
        s = []
        for key, value in values:
            s.append((self._orders.get(key, 99999), key, value))
        s.sort()
        return [(x, y) for z, x, y in s]

def process_value(value, encoding=None):
    length = len(value)
    t = value
    i = 0
    r = []
    buf = []
    listflag = False
    while i < length:
        if t[i] == '"': #string quote
            buf.append(t[i])
            i += 1
            while t[i] != '"' or (t[i] == '"' and t[i-1] == '\\'):
                buf.append(t[i])
                i += 1
            buf.append(t[i])
            i += 1
        elif t[i] == ',':
            r.append(''.join(buf))
            buf = []
            i += 1
            listflag = True
        elif t[i] == 'u':
            buf.append(t[i])
            i += 1
        else:
            buf.append(t[i])
            i += 1
            while i < length and t[i] != ',':
                buf.append(t[i])
                i += 1
    if buf:
        r.append(''.join(buf))
    result = []
    for i in r:
        if i.isdigit():
            result.append(int(i))
        elif i and i.startswith('u"'):
            result.append(unicode(unescstr(i[1:]), encoding))
        else:
            b = r_float.match(i)
            if b:
                result.append(float(b.group()))
            else:
                result.append(unescstr(i))

    if listflag:
        return result
    elif result:
        return result[0]
    else:
        return ''
    
unescapechars = {'"':'"', 't':'\t', 'r':'\r', 'n':'\n', '\\':'\\'}
def unescstr(value):
    if value.startswith('"') and value.endswith('"'):
        s = []
        i = 1
        end = len(value) - 1
        while i < end:
            if value[i] == '\\' and unescapechars.has_key(value[i+1]):
                s.append(unescapechars[value[i+1]])
                i += 2
            else:
                s.append(value[i])
                i += 1
        value = ''.join(s)
    return value

escapechars = {'"':r'\"', '\t':r'\t', '\r':r'\r', '\n':r'\n', '\\':r'\\'}
def escstr(value):
    s = []
    for c in value:
        if escapechars.has_key(c):
            s.append(escapechars[c])
        else:
            s.append(c)
    return ''.join(s)
    
def getdefaultencoding(encoding):
    if not encoding:
        encoding = locale.getdefaultlocale()[1]
    if not encoding:
        encoding = sys.getfilesystemencoding()
    if not encoding:
        encoding = 'utf-8'
    return encoding

def uni_prt(a, encoding=None):
    s = []
    if isinstance(a, (list, tuple)):
        for i, k in enumerate(a):
            s.append(uni_prt(k, encoding))
            s.append(',')
    elif isinstance(a, str):
        t = escstr(a)
        if ' ' in t or ',' in t or t.isdigit():
            s.append('"%s"' % t)
        else:
            s.append("%s" % t)
    elif isinstance(a, unicode):
        t = escstr(a)
        s.append('u"%s"' % t.encode(encoding))
    else:
        s.append(str(a))
    return ''.join(s)

if __name__ == '__main__':
    d = DictIni('t.ini')
    d._comment = 'Test\nTest2'
    d.a = 'b'
    d.b = 1
    d['b'] = 3
    d.c.d = (1,2,'b asf aaa')
    d['s']['t'] = u'中国'
    d['s'].a = 1
    d['m/m'] = 'testing'
    d.t.m.p = '3'
    d.save()
    
#    p = DictIni()
#    p.x = 1
#    
#    d.q = [p, p]
#
#    x = d.dict()
#    print type(x), x
#    print type(x['t']['m']), x['t']['m']
#    d.setfilename('test.ini')
#    d.t.m.comment('p', 'PTesting')
#    print d.getfilename()
#    print '---------------------'
#
#    d.save('test.ini')
#
#    a = process_value('1,abc,"aa cc",,"  ,\\"sdf",u"aaa"', 'ascii')
#    print a
#    print uni_prt(a, 'utf-8')
#
#    t = DictIni(inifile='test.ini')
#    print t
#
#    t.setfilename('test1.ini')
#    t.save()

#    t.setfilename('test2.ini')
#    t.save()