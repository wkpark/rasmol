/* prepdoc.c
 * Document Preparation System
 * Roger Sayle, January 1994
 * Version 1.0
 */

#include <string.h>
#include <stdio.h>
#include <ctype.h>

#ifndef True
#define True  1
#define False 0
#endif

#define RTFForm     0x00
#define LaTeXForm   0x01
#define HTMLForm    0x02
#define TextForm    0x03
#define HelpForm    0x04
#define ManForm     0x05



static char Buffer[82];
static FILE *InFile;
static int Format;


static int ReadLine()
{
    register char *ptr;
    register int len;
    register int ch;

    if( feof(InFile) )
    {   *Buffer = 0;
        return( False );
    }

    ptr = Buffer;
    do {
        ch = getc(InFile);
        if( (ch=='\n') || (ch=='\r') )
        {   while( ptr != Buffer )
                if( *(ptr-1)!=' ' )
                {   *ptr = 0;
                    return( True );
                } else ptr--;
        } else if( ch==EOF )
        {   *ptr = 0;
            return( True );
        } else *ptr++ = ch;
    } while( ptr < Buffer+80 );
    *ptr = 0;

    /* skip to the end of the line! */
    do { ch = getc(InFile);
    } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );
    return( True );
}


static int TextFlag;
static int TextCol;

static void DisplayLowerCase( ptr )
    char *ptr;
{
    register char ch;

    while( (ch = *ptr++) )
        if( isupper(ch) )
        {   putchar(tolower(ch));
        } else putchar(ch);
}

static void DisplayOnlyLowerCase( ptr )
    char *ptr;
{
    register char ch;
    
    while( (ch = *ptr++) )
        if( isupper(ch) )
        {   putchar(tolower(ch));
        } else if( ch != ' ' )
            putchar(ch);
}
    

static void DisplayText( src )
    char *src;
{
    register char *ptr;
    register char *dst;

    while( *src )
    {   dst = src;
        while( *dst && *dst!=' ' )
            dst++;

        if( TextCol+(dst-src) > 76 )
        {   putchar('\n');
            TextCol = 0;
        }

        while( src != dst )
        {   putchar( *src++ );
            TextCol++;
        }

        while( *src==' ' )
            src++;
        putchar(' ');
        TextCol++;
    }
}

static void ProcessCommand()
{
    static char buffer[80];
    register char *ptr;
    register int i,len;

    ptr = Buffer+2;
    switch( Buffer[1] )
    {   case('R'):  if( TextCol )
                        putchar('\n');
                    printf("%s\n",ptr);
                    TextFlag = True;
                    TextCol = 0;
                    break;

        case('T'):  if( Format == RTFForm )
                    {   printf("%s {}\n",ptr);
                        TextCol = 0;
                    } else if( (Format!=TextForm) && (Format!=HelpForm) )
                    {   printf("%s\n",ptr);
                        TextCol = 0;
                    } else DisplayText(ptr);
                    TextFlag = True;
                    break;

        case('P'):  if( TextFlag )
                    {   if( Format==HTMLForm )
                        {   fputs("<p>\n",stdout);
                        } else if( Format==RTFForm )
                        {   fputs("\\par\\par\n",stdout);
                        } else if( (Format==TextForm) || (Format==HelpForm) )
                        {   fputs("\n\n",stdout);
                        } else putchar('\n');
                        TextFlag = False;
                        TextCol = 0;
                    }
                    break;

        case('B'):  if( TextFlag )
                    {   if( Format==HTMLForm )
                        {   fputs("<p><hr>\n",stdout);
                        } else if( Format==RTFForm )
                        {   fputs("\\par\\page\\par\n",stdout);
                            fputs("+{\\footnote doc}\n",stdout);
                        } else if( (Format==TextForm) || (Format==HelpForm) )
                        {   fputs("\n\n",stdout);
                        } else putchar('\n');
                        TextFlag = False;
                        TextCol = 0;
                    }
                    break;

        case('S'):  if( Format==TextForm )
                    {   len = strlen(ptr);
                        printf("%s\n",ptr);
                        for( i=0; i<len; i++ )
                            putchar('-');
                        fputs("\n\n",stdout);
                    } else if( Format==HTMLForm )
                    {   fputs("<a name=\"",stdout);
                        DisplayOnlyLowerCase(ptr);
                        printf("\"><h3>%s</h3></a><p>\n",ptr);
                    } else if( Format==RTFForm )
                    {   fputs("#{\\footnote ",stdout);
                        DisplayOnlyLowerCase(ptr);
                        printf("}\n${\\footnote %s}\nK{\\footnote ",ptr);
                        DisplayLowerCase(ptr);
                        printf("}\n{\\b %s}\\par\\par\n",ptr);
                    } else if( Format==HelpForm )
                    {   putchar('?');
                        DisplayLowerCase(ptr);
                        printf("\n%s\n",ptr);
                    } else if( Format==ManForm )
                    {   printf(".TP\n.B %s\n",ptr);
                    }
                    TextCol = 0;
                    break;

        case('X'):  while( *ptr!=' ' )
                        ptr++;
                    *ptr++ = 0;
                    
                    if( Format == RTFForm )
                    {   printf("{\\uldb %s}{\\v %s} {}\n",ptr,Buffer+2);
                        break;
                    } else if( Format == HTMLForm )
                    {   printf("<a href=\"#%s\"><tt><b>%s</b></tt></a>\n",
                               Buffer+2,ptr);
                        break;
                    }
                    
        case('C'):  if( Format == RTFForm )
                    {   if( *ptr=='"' )
                        {   printf("\"{\\f2\\b %s}\" {}\n",ptr+1);
                        } else printf("{\\f2\\b %s} {}\n",ptr);
                        TextCol = 0;
                    } else if( Format == HTMLForm )
                    {   if( *ptr=='"' )
                        {   printf("\"<tt><b>%s</b></tt>\"\n",ptr+1);
                        } else printf("<tt><b>%s</b></tt>\n",ptr);
                    } else if( Format == ManForm )
                    {   if( *ptr=='"') ptr++;
                        printf(".B %s\n",ptr);
                    } else if( (Format!=TextForm) && (Format!=HelpForm) )
                    {   if( *ptr=='*' )
                        {   printf("\"%s\"\n",ptr);
                        } else printf("`%s'\n",ptr);
                        TextCol = 0;
                    } else /* DisplayText! */
                    {   if( *ptr=='"' )
                        {   sprintf(buffer,"\"%s\"",ptr+1);
                        } else sprintf(buffer,"`%s'",ptr);
                        DisplayText(buffer);
                    }
                    TextFlag = True;
                    break;
    }
}


int main( argc, argv )
    int argc;  char *argv[];
{
    register char *fname;
    register char *ptr;
    register int flag;

    fputs("Document Preparation System\n",stderr);
    fputs("Roger Sayle, January 1994\n",stderr);
    fputs("Version 1.0\n\n",stderr);

    Format = TextForm;

    if( argc==2 ) 
    {   fname = argv[1];
    } else if( argc==3 )
    {   fname = argv[2];
        ptr = argv[1];

        if( *ptr=='-' )
            ptr++;

        if( !strcmp(ptr,"latex") )
        {   Format = LaTeXForm;
        } else if( !strcmp(ptr,"html") )
        {   Format = HTMLForm;
        } else if( !strcmp(ptr,"help") )
        {   Format = HelpForm;
        } else if( !strcmp(ptr,"rtf") || !strcmp(ptr,"mshelp") )
        {   Format = RTFForm;
        } else if( !strcmp(ptr,"text") || !strcmp(ptr,"ascii") )
        {   Format = TextForm;
        } else if( !strcmp(ptr,"man") || !strcmp(ptr,"troff") )
        {   Format = ManForm;
        } else
        {   fputs("Formats:  -latex  LaTeX .tex file\n",stderr);
            fputs("          -troff  UNIX man(1) pages\n",stderr);
            fputs("          -html   HyperText metalanguage\n",stderr);
            fputs("          -help   RasMol on-line help file\n",stderr);
            fputs("          -rtf    Microsoft Help (Rich Text)\n",stderr);
            fputs("          -text   Standard ASCII text\n\n",stderr);
            exit(1);
        }

    } else /* DisplayUsage */
    {   fputs("Usage: prepdoc [format] <filename>\n",stderr);
        exit(1);
    }

    if( !(InFile=fopen(fname,"r")) )
    {   fputs("Error: Unable to open input file!\n",stderr);
        exit(1);
    }

    TextFlag = False;
    TextCol = 0;

    while( !feof(InFile) )
    {   ReadLine();
        switch( *Buffer )
        {   case('V'):  flag = (Format==LaTeXForm) ||
                               (Format==HTMLForm)  ||
                               (Format==RTFForm)   ||
                               (Format==TextForm);   break;

            case('D'):  flag = (Format==TextForm) ||
                               (Format==HelpForm) ||
                               (Format==ManForm);    break;

            case('S'):  flag = (Format==HelpForm) ||
                               (Format==ManForm);    break;

            case('N'):  flag = (Format==TextForm) ||
                               (Format==HelpForm);   break;

            case('L'):  flag = (Format==LaTeXForm);  break;
            case('H'):  flag = (Format==HTMLForm);   break;
            case('P'):  flag = (Format==HelpForm);   break;
            case('T'):  flag = (Format==TextForm);   break;
            case('M'):  flag = (Format==ManForm);    break;
            case('R'):  flag = (Format==RTFForm);    break;
            case('A'):  flag = True;                 break;
            default:    flag = False;
        }

        if( flag )
            ProcessCommand();
    }
    fclose(InFile);
    exit(0);
}

