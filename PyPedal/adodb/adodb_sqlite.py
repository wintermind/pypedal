########################################################################
# Vers 2.01 5 May 2006, (c)2005 GlennWashburn (crass#berlios.de) and John Lim
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
#
#  Requires http://initd.org/tracker/pysqlite
# Currently, the host, user, and password connection parameters are ignored and
# Affected_Rows() does not return the correct value
########################################################################
import adodb

try:
    import sqlite3 as sqlite
except:
    import pysqlite2.dbapi2 as sqlite
    
try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0


# Thread Safety= 2  connections
# Param Style  = pyformat "%(name)s"

class adodb_sqlite(adodb.ADOConnection):
    databaseType = 'sqlite'
    dataProvider = 'sqlite'
    
    sysDate = "date('now')"
    sysTimeStamp = "(date('now') || ' ' || time('now'))"

    metaColSQL = """PRAGMA table_info(%s)"""
    
    def __init__(self):
        pass

    def Module(self):
        return sqlite
    
    def _connect(self,host=None,user=None,password=None,database=None):
        # sqlite doesn't use host, user, or password
        dsn = database
        self._conn = sqlite.connect(dsn, detect_types=sqlite.PARSE_DECLTYPES)
        self._autocommit_conn(1)

    def _newcursor(self,rs):
        return cursor_sqlite(rs,self)
    
    def _autocommit_conn(self, level=1):
        if level:
            self._prev_iso_level = self._conn.isolation_level;
            self._conn.isolation_level = None;
        else:
            self._conn.isolation_level = self._prev_iso_level;

    def SelectLimit(self,sql,limit,offset=-1,params=None):
        if (offset >= 0): offset = " OFFSET "+str(offset)
        else: offset = ""
        return self.Execute(sql+" LIMIT "+str(limit)+offset,params)

    def BeginTrans(self):
        if self._autocommit:
            self._autocommit = False
        self._autocommit_conn(0)

    def RollbackTrans(self):
        self._conn.rollback()
        self._autocommit = True
        self._autocommit_conn(1)
        
    def CommitTrans(self):
        self._conn.commit()
        self._autocommit = True
        self._autocommit_conn(1)

    def UpdateBlob(self,table,field,blob,where,blobtype='BLOB'):
        self.Execute("update %s set %s = ? WHERE %s" % (table,field,where), (sqlite.Binary(blob),))
        
class cursor_sqlite(adodb.ADOCursor):
    def __init__(self,rs,conn):
        adodb.ADOCursor.__init__(self,rs,conn)
        self._insertid = rs.lastrowid

    
if __name__ == '__main__':
    db = adodb_sqlite()
    db.Connect('localhost','tester','test','test.sqlite')
    
    #~ db.Execute('create table adoxyz (id int, firstname text, lastname text, created timestamp default "2005-11-02");')
    #~ db.Execute("insert into adoxyz values (0, 'test1', 'ltest1', date(\"now\"));")
    #~ db.Execute("insert into adoxyz values (1, 'test2', 'ltest2', date(\"now\"));")
    #~ db.Execute("insert into adoxyz values (2, 'test3', 'ltest3', date(\"now\"));")
    #~ db.Execute("insert into adoxyz values (3, 'test4', 'ltest4', date(\"now\"));")
    
    #~ cur = db.Execute("select * from sqlite_master;")
    #~ print [r for r in cur];
    #~ cur = db.Execute("pragma table_info(adoxyz);")
    #~ print [r for r in cur];

    adodb.Test(db)
    
    #~ db.Execute('create table photos (id int, photo blob);')
    #~ db.Execute("insert into photos values (1, NULL);")
    adodb.Test_Blob(db)
