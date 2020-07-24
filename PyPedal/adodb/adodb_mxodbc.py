########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
#######################################################################

import adodb,adodb_odbc,sys

if sys.platform.find('win32') >= 0:
    import mx.ODBC.Windows
    mxODBC = mx.ODBC.Windows
else:
    import mx.ODBC.iODBC
    mxODBC = mx.ODBC.iODBC
    
try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0

class adodb_mxodbc(adodb.ADOConnection):
    databaseType = 'mxodbc'
    dataProvider = 'mxodbc'
    hasRowCount = False

    def __init__(self):
        pass

    def Module(self):
        global mxODBC
        
        return mxODBC
        #host=host1 user=user1 password=secret port=4341
    
    def _connect(self,host=None,user=None,password=None,database=None):
        global mxODBC
        
        if user == None and password == None and database == None:
            dsn = host
        else:
            dsn = 'dsn='+self.addq(host)+';'
            if (user != None): dsn += ' uid='+self.addq(user)+';'
            if (password != None): dsn += ' pwd='+self.addq(password)+';'
            if (database != None): dsn += ' database='+self.addq(database)+';'
      
        self._conn = mxODBC.DriverConnect(dsn)
        self._conn.setconnectoption(mxODBC.SQL.AUTOCOMMIT, mxODBC.SQL.AUTOCOMMIT_ON)

    def _newcursor(self,rs):
        return cursor_mxodbc(rs,self)
    
    def BeginTrans(self):
        global mxODBC
        if self._autocommit:
            self._autocommit = False
        self._conn.setconnectoption(mxODBC.SQL.AUTOCOMMIT, mxODBC.SQL.AUTOCOMMIT_OFF)

    def RollbackTrans(self):
        global mxODBC
        self._conn.rollback()
        self._autocommit = True
        self._conn.setconnectoption(mxODBC.SQL.AUTOCOMMIT, mxODBC.SQL.AUTOCOMMIT_ON)
        
    def CommitTrans(self):
        global mxODBC
        self._conn.commit()
        self._autocommit = True
        self._conn.setconnectoption(mxODBC.SQL.AUTOCOMMIT, mxODBC.SQL.AUTOCOMMIT_ON)

    def MetaColumns(self, table):
        curs = self._conn.cursor()
        curs.columns(None, None, table) # columns('%', '%', table) ?
        rs = self._newcursor(curs)
        arr = []
        table = table.upper()
        while not rs.EOF:
            if rs.fields[2].upper() == table:
                arr.append((rs.fields[3],rs.fields[5],rs.fields[6]))
            rs.MoveNext()
        return arr
    
class cursor_mxodbc(adodb.ADOCursor):
    def __init__(self,rs,conn):
        adodb.ADOCursor.__init__(self,rs,conn)
        rs.arraysize=10
            
if __name__ == '__main__':
    db = adodb_mxodbc()
    #db.Connect("Driver=Microsoft Visual FoxPro Driver;UID=;PWD=;SourceDB=D:\\inetpub\\adodb\\adoxyz.DBC;SourceType=DBC;Exclusive=No;BackgroundFetch=No;Collate=Machine;Null=Yes;Deleted=Yes;")
    db.Connect("Driver={Microsoft Access Driver (*.mdb)};Dbq=d:\\inetpub\\adodb\\northwind.mdb;Uid=Admin;Pwd=;")
    adodb.Test(db)        