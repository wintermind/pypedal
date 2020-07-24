########################################################################
# Vers 2.02 21 May 2007, (c)2004-2007 John Lim (jlim#natsoft.com.my) All Rights Reserved
# Released under a BSD-style license. See LICENSE.txt.
# Download: http://adodb.sourceforge.net/#pydownload
########################################################################

import adodb,adodb_mxodbc,datetime
from adodb import ADOConnection,ADOCursor

try:
    True, False
except NameError:
    # Maintain compatibility with Python 2.2
    True, False = 1, 0
        
class adodb_mxoracle(adodb_mxodbc.adodb_mxodbc):
    databaseType = 'mxoracle'
    dataProvider = 'mxodbc'
    sysDate = 'trunc(SysDate)'
    sysTimeStamp = 'SysDate'
    replaceQuote = "''"
    
    def _newcursor(self,rs):
        return cursor_mxoracle(rs,self)

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
    
        
class cursor_mxoracle(adodb_mxodbc.cursor_mxodbc):
    def __init__(self,rs,conn):
        adodb_mxodbc.cursor_mxodbc.__init__(self,rs,conn)

if __name__ == '__main__':
    db = adodb_mxoracle()
    db.Connect('sherkhan','scott','natsoft') #"dsn=sherkhan;uid=scott;pwd=natsoft")
    adodb.Test(db)