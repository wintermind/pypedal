########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################
import adodb
from adodb import ADOConnection,ADOCursor
import cx_Oracle

# threadsafety=2 (connections)
# paramstyle=named (:name)
try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0
    
class adodb_oci8(adodb.ADOConnection):
    databaseType = 'oci8'
    dataProvider = 'oci8'
    replaceQuote = "''"
    metaColSQL = "select cname,coltype,width from col where tname='%s' order by colno"    

    sysDate = 'trunc(SYSDATE)'
    sysTimeStamp = 'SYSDATE'

    NLS_DATE_FORMAT = 'YYYY-MM-DD'  ## To include time, use 'RRRR-MM-DD HH24:MI:SS'
    
    def __init__(self):
        pass

    def Module(self):
        return cx_Oracle    
  
    def _connect(self,host=None,user=None,password=None,database=None):
        if user == None and password == None and database == None:
            self._conn = cx_Oracle.connect('','',host)
        else:
            if host == None: self._conn = cx_Oracle.connect(user,password)
            else: self._conn = cx_Oracle.connect(user,password,host)

        self._query("ALTER SESSION SET NLS_DATE_FORMAT='"+self.NLS_DATE_FORMAT+"'")
		        

    def _newcursor(self,rs):
        if self._autocommit: self._conn.commit()
        rs = cursor_oci8(rs,self)
        if rs._isselect: rs._rowcount = -1 # oci8 does not return recordcount
        return rs
    
    def BeginTrans(self):
        if self._autocommit:
            self._autocommit = False

    def RollbackTrans(self):
        self._conn.rollback()
        self._autocommit = True
        
    def CommitTrans(self):
        self._conn.commit()
        self._autocommit = True    

    # offset not supported
    def SelectLimit(self,sql,limit,offset=-1,params=None):
        if offset == -1:
            return self.Execute("select * from ("+sql+") where rownum <= "+str(limit), params)
        else:
            raise StandardError, "SelectLimit does not support offset: " + sql

    def DBDate(self,d):
        if d == None: return 'null'
        return "TO_DATE(%s,'%s')" % (d.strftime(self.fmtDate),self.NLS_DATE_FORMAT)

    def DBTimeStamp(self,d):
        if d == None: return 'null'
        return "TO_DATE(%s,'%s')" % (d.strftime(self.fmtTimeStamp),self.NLS_DATE_FORMAT)

    def UpdateBlob(self,table,field,blob,where,blobtype='BLOB'):
        cursor = self._conn.cursor()
        #self.debug = True
        if blobtype == 'BLOB':
            bt = cx_Oracle.BLOB
        else:
            bt = cx_Oracle.CLOB
        cursor.setinputsizes(adodblob = bt)
        sql = "update %s set %s=:%s WHERE %s" % (table,field,'adodblob',where)
        self._query(sql,{'adodblob': blob},_cursor = cursor)

    def GetAll(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        if self.getLOBs:
            rows = []
            while 1:
                row = c.fetchone()
                if row == None: break
                rows.append(self._fixblobs(c.description,row))
            return rows
        else:
            return c.fetchall()

    def GetRow(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        arr = c.fetchone()
        return self._fixblobs(c.description,arr)

    def GetOne(self,sql,params=None):
        c = self._query(sql,params)
        if c == None: return None
        arr = c.fetchone()
        if (arr == None): return None
        if (c.description[0][1] == cx_Oracle.BLOB or c.description[0][1] == cx_Oracle.CLOB):
            return arr[0].read()
        return arr[0]

    def _hasblobs(self,description):
        for fld in description:
            t = fld[1]
            if t == cx_Oracle.BLOB or t == cx_Oracle.CLOB:
                return True
        return False
    
    def _fixblobs(self,description,row):
        if not row: return row
        
        arr = []
        i = 0
        for fld in description:
            t = fld[1]
            if t == cx_Oracle.BLOB or t == cx_Oracle.CLOB:
                try:
                    lob_fld=row.read()
                except:
                    arr.append(None)
                else:
                    arr.append(lob_fld)
            else:
                arr.append(row[i])
            i += 1
        return arr

    def MetaType(self, dtype):
        dtype = dtype.upper()
        if dtype == 'DATE':
            return 'T'
        return ADOConnection.MetaType(self, dtype)
    
    def MetaColumns(self, table):
        sql = self.metaColSQL % table.upper()
        return self.GetAll(sql)
        
class cursor_oci8(adodb.ADOCursor):
    getLOBs = None
    
    def __init__(self,rs,conn):
        if conn.getLOBs and rs.description:
        #    print rs.description
            self.getLOBs = conn._hasblobs(rs.description)
        else:
            self.getLOBs = False   
        ADOCursor.__init__(self,rs,conn)  

    def MoveNext(self):
        self.fields = self._cursor.fetchone()
        self.EOF = (self.fields == None)
        if self.getLOBs and not self.EOF:
            self.fields = self._conn._fixblobs(self._cursor.description,self.fields)
        return self.EOF

    def FetchRow(self):
        row = self.fields
        self.fields = self._cursor.fetchone()
        self.EOF = (self.fields == None)
        if self.getLOBs and not self.EOF:
            self.fields = self._conn._fixblobs(self._cursor.description,self.fields)
        return row      
        

if __name__ == '__main__':
    db = adodb_oci8()
    db.Connect('','scott','natsoft')
    adodb.Test(db)
    print db.MetaType('date')
    print db.MetaType('VarChar')
    print db.MetaType('DECiMAL')
    #adodb.Test_Blob(db)