########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################
import adodb
import psycopg

try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0

# Thread Safety= 2  connections
# Param Style  = pyformat "%(name)s"

class adodb_postgres(adodb.ADOConnection):
    databaseType = 'postgres'
    dataProvider = 'postgres'
    
    sysDate = "CURRENT_DATE"
    sysTimeStamp = "CURRENT_TIMESTAMP"

    metaColSQL = """SELECT a.attname,t.typname,a.attlen,a.atttypmod,a.attnotnull,a.atthasdef,a.attnum 
		FROM pg_class c, pg_attribute a,pg_type t 
		WHERE relkind = 'r' AND (c.relname='%s' or c.relname = lower('%s')) and a.attname not like '....%%'
AND a.attnum > 0 AND a.atttypid = t.oid AND a.attrelid = c.oid ORDER BY a.attnum"""
    
    def __init__(self):
        pass

    def Module(self):
        return psycopg
    
    #host=host1 user=user1 password=secret port=4341    
    def _connect(self,host=None,user=None,password=None,database=None):
        if user == None and password == None and database == None:
            dsn = host
        else:
            dsn = 'host='+self.addq(host)
            if (user != None): dsn += ' user='+self.addq(user)
            if (password != None): dsn += ' password='+self.addq(password)
            if (database != None): dsn += ' dbname='+self.addq(database)
        
        self._conn = psycopg.connect(dsn)
        self._conn.autocommit(1)

    def _newcursor(self,rs):
        return cursor_postgres(rs,self)

    def SelectLimit(self,sql,limit,offset=-1,params=None):
        if (offset >= 0): offset = " OFFSET "+str(offset)
        else: offset = ""
        return self.Execute(sql+" LIMIT "+str(limit)+offset,params)

    def BeginTrans(self):
        if self._autocommit:
            self._autocommit = False
        self._conn.autocommit(0)

    def RollbackTrans(self):
        self._conn.rollback()
        self._autocommit = True
        self._conn.autocommit(1)
        
    def CommitTrans(self):
        self._conn.commit()
        self._autocommit = True
        self._conn.autocommit(1)

    def _blobencode(self,blob):
        blob = str(blob)
        #92=backslash, 0=null, 39=single-quote
        return blob.replace(chr(92),r'\\134').replace(chr(0),r'\\000').replace(chr(39),r'\\047')

    def UpdateBlob(self,table,field,blob,where,blobtype='BLOB'):
        if (blobtype == 'BLOB'):
            self.Execute("update %s set %s='%s' WHERE %s" % (table,field,self._blobencode(blob),where))
        else:
            self.Execute("update %s set %s='%s' WHERE %s" % (table,field,self.addq(blob),where))
            
    def MetaColumns(self, table):
        #print self.metaColSQL
        sql = self.metaColSQL % (table,table)
        return self.GetAll(sql)            
        
class cursor_postgres(adodb.ADOCursor):
    def __init__(self,rs,conn):
        adodb.ADOCursor.__init__(self,rs,conn)

    
if __name__ == '__main__':
    db = adodb_postgres()
    db.Connect('localhost','tester','test','test')
    adodb.Test(db)
    #adodb.Test_Blob(db)