var MAX_RESULT = 20;

function updateSearch(text)
{
    if (DB == null) return;
    
    s = '';
    count = 1;
        
    if (text.length >= 2)
    {
        reg = new RegExp("^" + text.toLowerCase(), "gi");
        for (i = 0; i < DB.length-1; ++i)
        {
            entry = DB[i];
            key = entry[0];
            
            if (key.match(reg))
            {
                s += entry[1];
                ++count;
                if (count >= MAX_RESULT) break;
            }
        }
    }
        
    result.innerHTML = s ;
}
