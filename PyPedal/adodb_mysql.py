########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################

import adodb,re
import MySQLdb

#Thread Safety= 1 module
#Param Style  = format
 
try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0
    
class adodb_mysql(adodb.ADOConnection):
    databaseType = 'mysql'
    dataProvider = 'mysql'
    metaColSQL = "SHOW COLUMNS FROM %s"
    
    sysDate = 'CURDATE()'
    sysTimeStamp = 'NOW()'
	
   
    def __init__(self):
        pass

    def Module(self):
        return MySQLdb
    
    def _connect(self,host=None,user=None,password=None,database=None):
        self._conn = MySQLdb.connect(host, user, password, database)        

    def _newcursor(self,rs):
        return cursor_mysql(rs,self)

    def BeginTrans(self):
        if self._autocommit:
            self._autocommit = False
        self.Execute('set autocommit=0')
        #self._conn.Execute('begin')
        
    def RollbackTrans(self):
        self.Execute('rollback')
        self.Execute('set autocommit=1')
        self._autocommit = True

    def SelectLimit(self,sql,limit,offset=-1,params=None):
        if (offset >= 0): offset = str(offset)+","
        else: offset = ""
        return self.Execute(sql+" LIMIT "+offset+str(limit),params)
    
    def CommitTrans(self):
        self.Execute('commit')
        self.Execute('set autocommit=1')
        self._autocommit = True

    def UpdateBlob(self,table,field,blob,where,blobtype='BLOB'):
        self.Execute("update %s set %s='%s' WHERE %s" % (table,field,self.addq(blob),where))
 
    def qstr(self,s):
        return "'%s'" % self._conn.escape_string(s)

    def MetaColumns(self, table):
        sql = self.metaColSQL % table
        rs = self.Execute(sql)
        arr = []
        reFloat = re.compile("^(.+)\\((\\d+),(\\d+)")
        reInt = re.compile("^(.+)\\((\\d+)")
        while not rs.EOF:
            typeF = rs.fields[1]
            if typeF.find('(')>=0:
                if typeF.find(',')>=0:
                    m = reFloat.search(typeF) 
                else:
                    m = reInt.search(typeF)
                if m:
                    gps = m.groups()
                    type = gps[0]
                    size = gps[1]
                else:
                    type = typeF
                    size = -1
            else:
                type = typeF
                size = -1
                
            arr.append((rs.fields[0],type,size))
            rs.MoveNext()
        return arr
            
class cursor_mysql(adodb.ADOCursor):
    def __init__(self,rs,conn):
        adodb.ADOCursor.__init__(self,rs,conn)
        #self._insertid = rs.insert_id()
        self._insertid = rs.lastrowid    
if __name__ == '__main__':
    db = adodb.NewADOConnection('mysql')
    db.Connect('localhost','root','','northwind')
    for r in db.Execute('select * from adoxyz'):
        print r
    adodb.Test(db)