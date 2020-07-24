########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################

import adodb
import odbc,time,re,datetime

#Thread Safety= ? undefined
#Param Style  = ? undefined
try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0
        
class adodb_odbc(adodb.ADOConnection):
    databaseType = 'odbc'
    databaseProvider = 'odbc'
    hasRowCount = False
    replaceQuote = "''"
        
    def __init__(self):
        pass

    def Module(self):
        return odbc
    
    #host=host1 user=user1 password=secret port=4341    
    def _connect(self,host=None,user=None,password=None,database=None):
        if user == None and password == None and database == None:
            dsn = host
        else:
            dsn = 'dsn='+self.addq(host)
            if (user != None): dsn += '; uid='+self.addq(user)
            if (password != None): dsn += '; pwd='+self.addq(password)
            if (database != None): dsn += '; database='+self.addq(database)
        
        self._conn = odbc.odbc(dsn)
        self._conn.setautocommit(1)

    def _newcursor(self,rs):
        return cursor_odbc(rs,self)

    
    hasTop = 'TOP'
    _retop = None
    
    def SelectLimit(self,sql,limit,offset=-1,params=None):
        if self._retop == None:
            self._retop = re.compile(r'(^\s*select\s+(distinctrow|distinct)?)',re.IGNORECASE)
            
        sql = self._retop.sub('\\1 '+self.hasTop+' '+str(limit)+' ',sql)
        print "selectlimit: ",sql
        return self.Execute(sql)

    def BeginTrans(self):
        if self._autocommit:
            self._autocommit = False
        self._conn.setautocommit(0)

    def RollbackTrans(self):
        self._conn.rollback()
        self._autocommit = True
        
    def CommitTrans(self):
        self._conn.commit()
        self._autocommit = True

    # (1997, 12, 19, 1, 51, 53)
    def Date(self,d):
        return datetime.datetime(long(d[0]),long(d[1]),long(d[2]))
    
    def TimeStamp(self,d):
        d = time.localtime(d)
        return datetime.datetime(long(d[0]),long(d[1]),long(d[2]),long(d[3]),long(d[4]),long(d[5]))
    
class cursor_odbc(adodb.ADOCursor):
    def __init__(self,rs,conn):
        adodb.ADOCursor.__init__(self,rs,conn,norowcount=True)
            
if __name__ == '__main__':
    db = adodb_odbc()
    db.Connect("Driver={Microsoft Access Driver (*.mdb)};Dbq=d:\\inetpub\\adodb\\northwind.mdb;Uid=Admin;Pwd=;")
    adodb.Test(db)