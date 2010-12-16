#!/usr/bin/awk -f

BEGIN {

    LINE_FC_INDEX = 1
    LINE_T_INDEX = 1
#     split(substr($lines, 2, length($lines)), LINES_LIST)
    split(substr(lines_fc,2,length(lines_fc)), LINES_FC_LIST, "_")
    split(substr(lines_t,2,length(lines_t)), LINES_T_LIST, "_")
    NEWFC = 0 
    NEWT = 0

#     print "lines_fc " lines_fc | "cat 1>&2"
#     print "lines_fc_list[2] " LINES_FC_LIST[2] | "cat 1>&2"
#     print LINES_T_LIST[1] | "cat 1>&2"

}

NR == LINES_FC_LIST[LINE_FC_INDEX] {
    NEWFC = 1
}

NR == LINES_T_LIST[LINE_T_INDEX] {
    NEWT = 1
}


{
    if (NEWT == 1) {
        print $0 " //IOREV _todo_"
        LINE_T_INDEX = LINE_T_INDEX + 1
        NEWT = 0
    }
    else
        print $0
}



$0 ~ "{" {
#     print "reached braces"
    if (NEWFC == 1) {
        NEWFC = 0
        LINE_FC_INDEX = LINE_FC_INDEX + 1
        print "//IOREV _todo_"
    }
}
  