########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################

__author__ = "John Lim (jlim#natsoft.com)"
__credits__ = "(c) 2004-2007 John Lim"

import exceptions,sys,re
from datetime import datetime

try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0

MapTypes = {
            'VARCHAR' : 'C',
            'VARCHAR2' : 'C',
            'CHAR' : 'C',
            'C' : 'C',
            'STRING' : 'C',
            'NCHAR' : 'C',
            'NVARCHAR' : 'C',
            'VARYING' : 'C',
            'BPCHAR' : 'C',
            'CHARACTER' : 'C',
            'INTERVAL' : 'C',  # Postgres
            ##
            'LONGCHAR' : 'X',
            'TEXT' : 'X',
            'NTEXT' : 'X',
            'M' : 'X',
            'X' : 'X',
            'CLOB' : 'X',
            'NCLOB' : 'X',
            'LVARCHAR' : 'X',
            ##
            'BLOB' : 'B',
            'IMAGE' : 'B',
            'BINARY' : 'B',
            'VARBINARY' : 'B',
            'LONGBINARY' : 'B',
            'B' : 'B',
            ##
            'YEAR' : 'D', # mysql
            'DATE' : 'D',
            'D' : 'D',
            ##
            'TIME' : 'T',
            'TIMESTAMP' : 'T',
            'DATETIME' : 'T',
            'TIMESTAMPTZ' : 'T',
            'T' : 'T',
            ##
            'BOOL' : 'L',
            'BOOLEAN' : 'L', 
            'BIT' : 'L',
            'L' : 'L',
            ##
            'COUNTER' : 'R',
            'R' : 'R',
            'SERIAL' : 'R', # ifx
            'INT IDENTITY' : 'R',
            ##
            'INT' : 'I',
            'INTEGER' : 'I',
            'INTEGER UNSIGNED' : 'I',
            'SHORT' : 'I',
            'TINYINT' : 'I',
            'SMALLINT' : 'I',
            'I' : 'I',
            ##
            'LONG' : 'N', # interbase is numeric, oci8 is blob
            'BIGINT' : 'N', # this is bigger than PHP 32-bit integers
            'DECIMAL' : 'N',
            'DEC' : 'N',
            'REAL' : 'N',
            'DOUBLE' : 'N',
            'DOUBLE PRECISION' : 'N',
            'SMALLFLOAT' : 'N',
            'FLOAT' : 'N',
            'NUMBER' : 'N',
            'NUM' : 'N',
            'NUMERIC' : 'N',
            'MONEY' : 'N',
            ## informix 9.2
            'SQLINT' : 'I', 
            'SQLSERIAL' : 'I', 
            'SQLSMINT' : 'I', 
            'SQLSMFLOAT' : 'N', 
            'SQLFLOAT' : 'N', 
            'SQLMONEY' : 'N', 
            'SQLDECIMAL' : 'N', 
            'SQLDATE' : 'D', 
            'SQLVCHAR' : 'C', 
            'SQLCHAR' : 'C', 
            'SQLDTIME' : 'T', 
            'SQLINTERVAL' : 'N', 
            'SQLBYTES' : 'B', 
            'SQLTEXT' : 'X'}

class adodb_iter:
    cursor = None
    
    def __iter__(self):
        return self

    def next(self):
        if self.cursor.EOF: raise StopIteration
        ret = self.cursor.fields
        self.cursor.MoveNext()
        return ret

    

def NewADOConnection(modulename):

    if modulename.find(':') >= 0:
        # handle connection string of the form driver://user:pwd@server/database
        # where user, pwd, database are all optional
        match=re.match('([^:]*)://(.*)/([^\?]*)\?*(.*)', modulename)
        if match:
            gps = match.groups()
            server = ''
            user = ''
            pwd = ''
            db = ''

            if len(gps) >= 1:
                modulename = gps[0]

            if len(gps) >= 2:
                mid = gps[1]
                if mid.find('@') >= 0:
                    upwd,server = mid.split('@')
                    if mid.find(':') >= 0:
                        user,pwd = upwd.split(':')
                else:
                    if mid.find(':') >= 0:
                        user,pwd = mid.split(':')
                    else:
                        server = mid

            if len(gps) >= 3:
                db = gps[2]
                
            #print server, user, pwd, db
            conn = NewADOConnection(modulename)
            conn.Connect(server,user,pwd,db)
            return conn
        
    if modulename == 'oracle':
        modulename = 'oci8'
        
    try:
        modulename = 'adodb_'+modulename
        module = __import__(modulename,globals(), None, [modulename])
    except ImportError:
        return None
    klass = vars(module)[modulename]
    return klass()

def ADONewConnection(modulename):
    return NewADOConnection(modulename)

class ADOConnection:
    databaseType = None
    dataProvider = 'native'
    host = None
    user = None
    password = None
    database = None
    replaceQuote = "\\'"
    useExceptions = True
    debug = None
    getLOBs = True
    hasRowCount = True
    metaColSQL = 'Invalid'
    fmtDate = '%Y-%m-%d'
    fmtTimeStamp = '%Y-%m-%d %H:%M:%S'
    
    _errormsg = ''
    _errno = 0
    _conn = None
    _autocommit = True
    _connected = True
    
    def __init__(self):
        pass            
    
    def Connect(self,host=None,user=None,password=None,database=None):
        self.database = database
        self.host = host
        self.user = user
        self.password = password
        self._connect(host,user,password,database)
        return bool(self._conn)

    def IsConnected(self):
        return bool(self._conn)
    
    def DriverInfo(self):
        try:
            m = self.Module()
            print "Driver        =",self.databaseType
            print "API Level     =",m.apilevel
            print "Param Style   =",m.paramstyle
            print "Thread Safety =",m.threadsafety," (0=none, 1=module, 2=connections, 3=cursors)"
            print "--------------"
        except:
            print "???????"
    
    def ErrorMsg(self):
        return self._errormsg

    def ErrorNo(self):
        return self._errno
    

    def qstr(self,s):
        if (self.replaceQuote == "\\'"): s = str(s).replace('\\','\\\\')
        return "'"+str(s).replace("'", self.replaceQuote)+"'"
    
    def quote(self,s):
        return "'"+str(s).replace("'", self.replaceQuote)+"'"

    def addq(self,s):
        if (self.replaceQuote == "\\'"): s = str(s).replace('\\','\\\\')
        return str(s).replace("'", self.replaceQuote)

    def Conn(self):
        return self._conn

    def _query(self,sql,params=None,_cursor=None):

        try:
            if _cursor == None: _cursor = self._conn.cursor()
            if self.debug:
                s = "(%s): %s" % (self.databaseType, sql)
                if type(self.debug) is str:
                    try:
                        f = file(self.debug,'w')
                        f.write(s+"\n")
                        f.close()
                    except:
                        pass
                elif long(self.debug) == 2:
                    print "<hr>"+s.replace('&','&amp;').replace('<','&lt;')+"<hr>"
                else:
                    print s;
            #if debug
                    
            if params == None: _cursor.execute(sql)
            else: _cursor.execute(sql,params)
            self._errormsg = ''
            self._errno = 0
            
        except StandardError, err:
            self._errormsg = str(err)
            self._errno = -1
            if self.useExceptions: # use default python error handling
                raise sys.exc_info()[0] ,str(err)+': '+sql                      
            _cursor = None

        return _cursor

    def SelectLimit(self,sql,limit,offset=-1,params=None):
        pass
    
    def Execute(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        rs = self._newcursor(c)
        return rs

    def UpdateBlob(self,table,field,blob,where,blobtype='BLOB'):
        raise StandardError, 'UpdateBlob not supported'

    def UpdateBlobFile(self,table,field,filepath,where,blobtype='BLOB'):
        f = file(filepath, 'rb')
        data = f.read()
        f.close()
        self.UpdateBlob(table,field,data,where,blobtype)

    def UpdateClob(self,table,field,blob,where):
        self.UpdateBlob(table,field,blob,where,'CLOB')

    def GetRows(self,sql,params=None):
        return self.GetAll(sql,params)
    
    def GetArray(self,sql,params=None):
        return self.GetAll(sql,params)
    
    def GetAll(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        all = c.fetchall()
        c.close()
        return all
    
    def GetRow(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        c.close()
        return c.fetchone()

    def GetRow(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        row = c.fetchone()
        c.close()
        return row

    def GetOne(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        arr = c.fetchone()
        c.close()
        if (arr == None): return None        
        return arr[0]

    def GetCol(self, sql, params=None):
        rs = self.Execute(sql,params)
        arr = []
        while not rs.EOF:
            arr.append(rs.fields[0])
            rs.MoveNext()
        rs.Close()
        return arr

    def GetAssoc(self, sql, params=None):
        rs = self.Execute(sql,params)
        dict = {}
        if rs.EOF:
            return None
        
        if len(rs.fields) == 2:
            while not rs.EOF:
                dict[rs.fields[0]] = rs.fields[1]
                rs.MoveNext()
        elif len(rs.fields)>2:
            while not rs.EOF:
                dict[rs.fields[0]] = rs.fields[1:]
                rs.MoveNext()
        else:
            while not rs.EOF:
                dict[rs.fields[0]] = None
                rs.MoveNext()
        
        return dict

    def GetDict(self, sql, params=None):
        return self.GetAssoc(sql,params)
    
    def BeginTrans(self):
        pass

    def CommitTrans(self):
        pass
    
    def RollbackTrans(self):
        pass
    
    def Close(self):
        try:
            if self._conn != None: self._conn.close()
        except:
            pass
        
        self._conn = None

    def DBDate(self,d):
        return "'%s'" % d.strftime(self.fmtDate)

    def DBTimeStamp(self,d):
        return "'%s'" % d.strftime(self.fmtTimeStamp)   

    _redate = None # compiled regex

    # convert iso format date to python timestamp    
    def Date(self,s):
        dates = str(s)
        if self._redate == None:
            self._redate = re.compile("^([0-9]{4})[-/\.]?([0-9]{1,2})[-/\.]?([0-9]{1,2})")
        match = self._redate.search(dates)
        if not match: return None
        
        year = long(match.group(1))
        month = long(match.group(2))
        day = long(match.group(3))
        return datetime(year, month, day)
    
    _rets = None # compiled regex

    # convert iso format timestamp to python timestamp    
    def TimeStamp(self,s):
        ts = str(s)
        if self._rets == None:
            self._rets = re.compile("^([0-9]{4})[-/\.]?([0-9]{1,2})[-/\.]?([0-9]{1,2})[ ,-]*([0-9]{1,2}):?([0-9]{1,2}):?([0-9\.]{1,5})")
        match = self._rets.search(ts)
        if not match: return self.Date(ts)
        
        year = long(match.group(1))
        month = long(match.group(2))
        day = long(match.group(3))
        hour = long(match.group(4))
        min = long(match.group(5))
        sec = long(float(match.group(6).strip())) # Python type-nazis
        return datetime(year, month, day, hour, min, sec)

    def MetaType(self, dbtype):
        global MapTypes
        
        dbtype = dbtype.upper()
        if MapTypes.has_key(dbtype):
            return MapTypes[dbtype]
      
        return 'N'

    def MetaColumns(self, table):
        sql = self.metaColSQL % table
        return self.GetAll(sql)

    
        
class ADOCursor:
    _cursor = None
    fields = None
    EOF = False
    _rowcount = 0
    _isselect = False
    _insertid = 0
    _conn = None
    
    def __init__(self,rs,conn,norowcount=False):
        self._cursor = rs
        self._conn = conn
        if norowcount: self._rowcount = -1
        else: self._rowcount = rs.rowcount
       
        try:
            self.MoveNext()
            self._isselect = True
        except:
            pass # not a select statement
        self.EOF = (self.fields == None)

    def __iter__(self):
        iter = adodb_iter()
        iter.cursor = self
        return iter
    
    def RecordCount(self):
        return self._rowcount
        
    def MoveNext(self):
        self.fields = self._cursor.fetchone()
        self.EOF = (self.fields == None)
        return self.EOF

    def FetchRow(self):
        row = self.fields
        self.fields = self._cursor.fetchone()
        self.EOF = (self.fields == None)
        return row

    # returns a tuple of the form (name, type_code,display_size, internal_size, precision, scale,null_ok)
    # note: databases could return name in upper or lower-case
    def FetchField(self,row):
        #print self._cursor.description
        if len(self._cursor.description) <= row: return None
        return self._cursor.description[row]
    
    def Affected_Rows(self):
        if self._rowcount >= 0: return self._rowcount
        return 0

    def Insert_ID(self):
        return self._insertid

    def Cursor(self):
        return self._cursor

    def GetRowAssoc(self,upper=1):
        d = {}
        i = 0
        desc = self._cursor.description
        if upper:
            for i in xrange(0,len(self.fields)):
                d[desc[i][0].upper()] = self.fields[i]
        elif not upper:
            for i in xrange(0,len(self.fields)):
                d[desc[i][0].lower()] = self.fields[i]
        else:
            for i in xrange(0,len(self.fields)):
                d[desc[i][0]] = self.fields[i]

        return d
    
    def Close(self):
        if self._cursor: self._cursor.close()
        self._cursor = None
        
#===========================================================
#                  UNIT TESTING
#===========================================================

def _Test_Eq(testid, correct, testval, errmsg=''):
    if correct == testval:
        print "Passed Test: "+testid
    else:
        print ""
        print "********* Failed Test: "+testid
        print "********************** "+str(errmsg)
        print "********************** expected="+str(correct)
        print "**********************   actual="+str(testval)

def Test_Blob(db):
    import os
    src = 'c:/lensserver.gif'
    dest = 'c:/testpy1.gif'
    try:
        os.unlink(dest)
    except:
        pass

    saved = db.debug
    saveb = db.getLOBs
    db.debug = True
    db.getLOBs = True
    
    db.UpdateBlobFile('photos','photo',src,'id=1')
    data = db.GetOne('select photo from photos where id=1')
    f = file(dest,'wb')
    f.write(data)
    f.close()

    rs = db.Execute('select * from photos')
    while not rs.EOF:
        print 'Fields=',rs.fields
        rs.MoveNext()
        
    print "======================="

    rows = db.GetAll('select * from photos where id<=1')
    print rows
    
    db.getLOBs = saveb
    db.debug = saved

def Test(db,debug=False):
    db.DriverInfo()

    if False:
        d = db.Date('2004-03-21')
        print '2004-03-21=',d

        d = db.TimeStamp('2004-03-22 12:50:51')
        print '2004-03-22 12:50:51=',d

        print "DBTimeStamp=", db.DBTimeStamp(d)
    
    db.useExceptions = True # use adodb error handling
    try:
        sql = 'select * from xadoxyz where 0 < id and id < 3'
        rs = db.Execute(sql)
        _Test_Eq('Bad SQL',None, rs, sql)
    except:
        print "And you should see an error message indicating bad table was defined: "
        print "err=",db.ErrorMsg()
    
    print "-----"
    rs = db.Execute('select * from ADOXYZ where 0 < id and id < 3 order by id')
    while not rs.EOF:
        print rs.fields
        rs.MoveNext()
        
    print "You should see 2 rows of data here:"
    rs = db.Execute('select * from adoxyz where 0 < id and id < 3 order by id')
    print "rows=",rs.RecordCount()
    while (not rs.EOF):
        print rs.GetRowAssoc()
        rs.MoveNext()

    print "-----"
    rs = db.Execute('select id,firstname from adoxyz where 0 < id and id < 3 order by id')
    _Test_Eq("Test FetchField",'FIRSTNAME',rs.FetchField(1)[0].upper())
    if (debug): print rs.FetchField(1)
    cnt = 0
    while 1:
        arr=rs.FetchRow()
        if arr == None: break
        cnt += 1
        _Test_Eq('Execute 2.0',cnt,arr[0])

    _Test_Eq('Execute 2.1',2,cnt)
    if rs.RecordCount() == -1: print "*** RecordCount not supported: -1"
    else: _Test_Eq('Execute 2.1 RecordCount',2,rs.RecordCount())  
    
    rs = db.Execute("delete from adoxyz where id=997")
    cnt = rs.Affected_Rows()
    _Test_Eq('Affected_Rows',1,cnt)
    
    ok = db.Execute("insert into adoxyz (id, firstname,lastname) values (997,'python','snake')")
    if not ok:  _Test_Eq('DELETE/INSERT','inserted row','failed insert')
       
    row = db.GetRow("select id,firstname from adoxyz where id=997");        
    _Test_Eq('GetRow',str(997)+' '+'python',str(int(row[0]))+' '+row[1].rstrip(),row)
    
    row = db.GetOne("select id,firstname from adoxyz where id=997");        
    _Test_Eq('GetOne',997,row)

    rs = db.SelectLimit("select id,firstname from adoxyz",3)
    cnt = 0
 
    try:
        for row in rs:
            cnt += 1
            #print rs.fields
        _Test_Eq('SelectLimit',3,cnt)
    except:
        print "Failed Iteration"
        print sys.exc_info()[1]
    
    d = db.GetOne('select created from adoxyz where id=1')
    d2 = db.TimeStamp(d)
    _Test_Eq('DBDate',str(d)[:19],str(d2))

    if (db.qstr("\\show'boat") != "'\\\\show\\'boat'" and db.qstr("\\show'boat") != "'\\show''boat'"):
        _Test_Eq('qstr',"qstr(\\show'boat)", db.qstr("\\show'boat"))
    else:
        _Test_Eq('qstr','1','1')

    try:
        db.debug=True
        print "Testing GetAssoc"
        arr = db.GetAssoc('select firstname,lastname from adoxyz')
        print arr
        print "Testing GetCol"
        arr = db.GetCol('select firstname from adoxyz')
        print arr
    except:
        print sys.exc_info()[1]
        
    try:
        print "MetaColumns:"
        rows = db.MetaColumns('adoxyz')
        print rows
    except:
        print "Failed MetaColumns"
        print sys.exc_info()[1]
    
    try:    
        db.BeginTrans()
        ok = db.Execute("insert into adoxyz (id, firstname,lastname) values (1997,'python','snake')")
        db.RollbackTrans()
        val = db.GetOne('select * from adoxyz where id=1997')
        _Test_Eq('Rollback Test',None,val)
    except:
        print "Failed Rollback Test"
        print sys.exc_info()[1]