########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################

import adodb,adodb_mxodbc,datetime

try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0
        
class adodb_access(adodb_mxodbc.adodb_mxodbc):
    databaseType = 'access'
    dataProvider = 'mxodbc'
    sysDate = "FORMAT(NOW,'yyyy-mm-dd')"
    sysTimeStamp = 'NOW'
    
    def _newcursor(self,rs):
        return cursor_access(rs,self)    
        
class cursor_access(adodb_mxodbc.cursor_mxodbc):
    def __init__(self,rs,conn):
        adodb_mxodbc.cursor_mxodbc.__init__(self,rs,conn)
        self._rowcount = rs.rowcount

if __name__ == '__main__':
    db = adodb_access()
    db.Connect("Driver={Microsoft Access Driver (*.mdb)};Dbq=d:\\inetpub\\adodb\\northwind.mdb;Uid=Admin;Pwd=;")
    adodb.Test(db)