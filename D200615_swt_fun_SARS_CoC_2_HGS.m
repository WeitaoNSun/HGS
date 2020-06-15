function D200615_swt_fun_SARS_CoC_2_HGS()
clear all;
close all hidden;
%set blast program folder;
dir_blast='D:\blast\blast-2.10.0+\';
%set human genome database folder;
dir_db='D:\blast\blast-2.10.0+\db\HumanG+T\';
%set output data folder;
paths='d:\mydata\';
addpath(strcat([dir_blast,'bin\'])); %add blast local excutive file folfer;
addpath(strcat([dir_blast,'db\']));  %add blast local library folder;
enum_data_id={'sars2_china_GISAID','sars2_usa_GISAID','sars2_Europe_GISAID'};
data_id=enum_data_id{1};
switch data_id
    case enum_data_id{1} %China
        filelst='D:\AAA_swt_portablework\mypapers\2020\2020.1.28-sars-spike\GISAID\AAAAA_GISAID_COVID19_CHINA_USA_ERUOPE\D200429-SARS-ncbi-149\fasta-no1a\sars-fsa.lst';
    case enum_data_id{2} %USA
        filelst='D:\AAA_swt_portablework\mypapers\2020\2020.1.28-sars-spike\GISAID\AAAAA_GISAID_COVID19_CHINA_USA_ERUOPE\SARS2_USA\fsa\sars2_usa_fsa.lst';
    case enum_data_id{3} %Europe
    otherwise
        return;
end
virus_tag='SARS-CoV-2';
readme={'------------ WeitaoSUN@TsinghuaUniv ---------',...
    '*.hgs: blastn marks, expection values, and HGS values for each segment/orf of virus strain.',...
    '*.txt: details of segment/orf matching information in database.'};

file_strain_CDS='d:\mydata\EPI_ISL_431180_norf10.fsa';
[~,filetag,ext]=fileparts(file_strain_CDS);
data_tag=strcat([virus_tag,'_',filetag]);
[accession]=swt_fun_dataset_HGS(file_strain_CDS,data_tag);

%-------------------------------------------------------------------------
    function [accession]=swt_fun_dataset_HGS(file_strain_CDS,data_tag)
        file_ncbi_virus_fasta_CDS=strcat([file_strain_CDS]);
        % file_ncbi_virus_fasta_CDS='D:\mydata\all_matlab_picture\D20200129_ncov\20200418\AY394850.2_ORF1ab.fsa';
        
        %------------read fasta file containing all virus strain ORF sequences--------;
        [Header, Sequence] = fastaread(file_ncbi_virus_fasta_CDS);%,'blockread', [1 5]);
        if(~iscell(Sequence))
            Header={Header};
            Sequence={Sequence};
            nseq=1;
        else
            nseq=length(Sequence); %count CDS sequence number of current virus;
        end
        %get Access number and other information from fasta header section;
        if(~iscell(Header))
            header_split=split(Header,'|');
        else
            header_split=cellfun(@(x) split(x,'|'), Header,'un',0);
        end
        %get access number;
        accession=cellfun(@(v)v(1),header_split);
        
        %         nstrain=length(unique(AnumAll));  %count accession numbers of each segment/orf of current virus;
        
        % header_split=cellfun(@(x) split(x,':'), Header,'un',0);
        % Anum=cellfun(@swt_fun_fasta_header_Anum,header_split,'un',0);
        % ORF=cellfun(@swt_fun_fasta_header_ORF,header_split,'un',0);
        % ORF_range=cellfun(@swt_fun_fasta_header_ORF_range,header_split,'un',0);
        %---------------------make Human G+T database ------------------------
        humanG38=strcat([dir_db,'ncbi-genomes-2020-04-15\GCF_000001405.38_GRCh38.p12_genomic.fna']);
        
        %-----------------------initialized HGS output files ---------------------------------------
        str_date=datestr(now,'yyyymmdd');
        file_tag=strcat([str_date,'_(',data_tag,')_seq',num2str(nseq)]);
        file_HGS=strcat([paths,file_tag,'.hgs']);
        fid=fopen(file_HGS,'w');
        % fprintf(fid,'AccessNum\tname\tloc\tstart\tend \tSUM(S)\tseqlen\tHGS\r\n');
        fprintf(fid,'max(S)\tsum(S)\tseqlen\tHGS\tquery\tsequence\r\n');
        fclose(fid);
        
        str_date=datestr(now,'yyyymmdd');
        file_log=strcat(paths,file_tag,'.log');
        
        %----------- load each sequence and calculate HGS --------------------
        is_seq_size_limit=0; seq_limit=3000;  %if sequence is too long, skip it. (option)
        format_id=6;
        str_blastargs=strcat([' -W 11  -T T -r 2 -q 3 -culling_limit 1 -m ',num2str(format_id)]);
        %   -a  Number of processors to use [Integer],     default = 1
        %see following help on blast command.
        %   -r  Reward for a nucleotide match (blastn only) [Integer],      default = 1
        %   -W  Word size, default if zero (blastn 11, megablast 28, all others 3) [Integer],     default = 0
        %   -K  Number of best hits from a region to keep (off by default, if
        %   used a value of 100 is recommended) [Integer],      default = 0
        %   -q  Penalty for a nucleotide mismatch (blastn only) [Integer]
        %     default = -3
        %   -r  Reward for a nucleotide match (blastn only) [Integer]
        %     default = 1
        %   -T  Produce HTML output [T/F],    default = F
        %         switch format_id
        %             case {1,2,3,4,5,6}
        %                 %---------- create HGS result file for current sequence -------------
        %                 %             fileout=strcat(['blastn_',data_tag,'_',Header{i},'_',db_title,'_len',num2str(seq_len),'_fmt',num2str(format_id),'.txt']);
        %                 %                     fileout=strcat(['blastn_',data_tag,'_',db_title,'_len',num2str(seq_len),'_fmt',num2str(format_id),'.txt']);
        %                 fileout=strcat(['blastn_',data_tag,'_',db_title,'_fmt',num2str(format_id),'.txt']);
        %                 fileout=strrep(fileout,'|','-');      fileout=strrep(fileout,':','-'); %remove odd chars in file name;
        %             otherwise
        %                 %---------- create HGS result file for current sequence -------------
        %                 %             fileout=strcat(['blastn_',data_tag,'_',Header{i},'_',db_title,'_len',num2str(seq_len),'_fmt',num2str(format_id),'.txt']);
        %                 %                     fileout=strcat(['blastn_',data_tag,'_',db_title,'_len',num2str(seq_len),'_fmt',num2str(format_id),'.txt']);
        %                 fileout=strcat(['blastn_',data_tag,'_',db_title,'_fmt',num2str(format_id),'.txt']);
        %                 fileout=strrep(fileout,'|','-');      fileout=strrep(fileout,':','-'); %remove odd chars in file name;
        %         end
        fileout=strcat([paths,file_tag,'_fmt',num2str(format_id),'.txt']);
        
        fileout=strrep(fileout,' ','-'); %remove blanks in filename, if there are blank, command line will fail.
        fileout=strrep(fileout,'/','-'); %remove blanks in filename, if there are blank, command line will fail.
        fileout=strrep(fileout,'>','');
        fid=fopen(fileout,'w');
        fprintf(fid,'query\tTarget\tPercetIdent\tAlignmentLen\tMismatchNo\tgapNo\tQueryStart\tQueryEnd\tTargetStart\tTargetEnd\tEValue\tBitScore\r\n');
        fclose(fid);
        %          Field	 	Description
        % 1	 	Query label.
        % 2	 	Target (database sequence or cluster centroid) label.
        % 3	 	Percent identity.
        % 4	 	Alignment length.
        % 5	 	Number of mismatches.
        % 6	 	Number of gap opens.
        % 7	 	Start position in query. Query coordinates start with 1 at the first base in the sequence as it appears in the input file. For translated searches (nucleotide queries, protein targets), query start<end for +ve frame and start>end for -ve frame.
        % 8	 	End position in query.
        % 9	 	Start position in target. Target coordinates start with 1 at the first base in sequence as it appears in the database. For untranslated nucleotide searches, target start<end for plus strand, start>end for a reverse-complement alignment.
        % 10	 	End position in target.
        % 11	 	E-value calculated using Karlin-Altschul statistics.
        % 12	 	Bit score calculated using Karlin-Altschul statistics.
        
        for i=1:nseq
            %---- get access number and length of current sequence -------
            headeri=split(Header{i},'|');
            Anumi=headeri{1};
            ORFi=headeri{2};
            seq_len=length(Sequence{i});
            
            if(is_seq_size_limit)
                if(seq_len>seq_limit)
                    continue;
                end
            end
            %some times there are strang chars in sequence, such as "NNNNNN" in MT259279.1	ORF6. we
            %skip such sequences;
            odd_nucleotide='NNN';
            idx_strage_seq= (strfind(Sequence{i},odd_nucleotide));
            if(~isempty(idx_strage_seq))
                % there are straing sequence character;
                a=0;
                fid=fopen(file_log,'a');
                fprintf(fid,'Warning:\tthere are %s in sequence %s:\t%s\t%s\r\n',odd_nucleotide,Anumi,Header{i},Sequence{i});
                fclose(fid);
                continue;
            else
                a=0;
                %  continue;
            end
            
            %save current sequence to a temp fasta file;
            filefast=strcat([paths,'temp.fsa']);
            fid=fopen(filefast,'w');
            fclose(fid);
            filefastout=strcat([paths,'temp.txt']);
            fido=fopen(filefastout,'w');
            fclose(fido);
            
            fastawrite(filefast,strrep(Header{i},' ',''),Sequence{i});
            
            db_used=humanG38;
            db_title='humanG38';
            if(~isfile(db_used))
                display(strcat(['Data base ',db_used,' does not exists.']));
                return;
            end
            %---------- run HGS calculation ----------------------------
            %     You are welcome. To search against multiple BLAST databases, just concatenate the database names separated by space
            % -db "nr swissprot"
            %     the default ¡®expect threshold¡¯ value E=10 by BLAST system is used, and the seed sequence length is 11. Match/Mismatch score is set as 2/-3. Gap costs are chosen as: existence=5 and extension=2.
            %     result=blastlocal('inputquery', filefast, 'Program','blastn','database',db_used,...
            %         'tofile',fileout, 'Filter',true, 'GapOpen',-5, 'GapExtend',-2,...
            %        'blastargs',str_blastargs);
            
            %------------------ Main Part -----------------------------
            out=swt_fun_blastn(filefast,db_used,filefastout,format_id);
            
            %--------------------------------------------------------------
            %save HGS data;
            [maxS,TotalS,HGS,Accession]=swt_fun_blastn_HGS_to_file(seq_len,filefastout,file_HGS,Anumi,Sequence{i});
            system(['type ',filefastout,' >>',fileout]);
            swt_fun_write_text(fileout,'');
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxS,TotalS,HGS,Accession]=swt_fun_blastn_HGS_to_file(len,fileout,fileHGS,Anum,sequence)
maxS=[];
TotalS=[];
HGS=[];
Accession=cell({});
if(~isfile(fileout))
    return
end
data=importdata(fileout);
if(isempty(data))
    return;
    fid=fopen(fileHGS,'a');
    fprintf(fid,'%s\t%.2f\t%.2f\t%d\t%.4f\r\n',Anum,NaN,NaN,NaN,NaN);
    fclose(fid);
    return;
end
Access=data.textdata(:,2);
AccessUniq=unique(Access);
nuniq=length(AccessUniq);
for i=1:nuniq
    idx=find(ismember(Access,AccessUniq(i)));
    query(i)=data.textdata(idx(1));
    Accession(i)=data.textdata(idx(1),2);
    Si=data.data(idx,end);
    maxS(i,1)=max(Si);
    TotalS(i,1)=sum(Si);
end
HGS=sum(TotalS)/(2*len);

fid=fopen(fileHGS,'a');
fprintf(fid,'%.2f\t%.2f\t%d\t%.4f\t%s\t%s\t\r\n',max(maxS),sum(TotalS),len,HGS,data.textdata{1},sequence);
fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function swt_fun_blastn_result(result,seq_len,fileout,is_verbose)
if(isempty(result))
    return;
end
query=result.Query;
nHits=length(result.Hits);
for i=1:nHits
    nHSP=length(result.Hits(i).HSPs);
    Accession{i,1}=result.Hits(i).Name;
    scores=[result.Hits(i).HSPs(:).Score];
    Evalue=[result.Hits(i).HSPs(:).Expect];
    Ident=[result.Hits(i).HSPs.Identities];
    PerIdent=[Ident.Percent];
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    MaxScore(i,1)=max(scores);
    TotalScore(i,1)=sum(scores);
    PerIdent_show(i,1)=max(PerIdent);  %only show the largest(best) PerIdent;  The larger this value, the better;
    Evalue_show(i,1)=min(Evalue); %only show the smallest(best) Evalue, the smaller this value, the better;
end
HGS=sum(MaxScore)/(2*seq_len);
header=split(result.Query,'|');

%wirte HGS data file;
fid=fopen(fileout,'a');
if(is_verbose) %------------------------all information format --------------------
    for i=1:nHits
        fprintf(fid,'%s\t%s\t%.1f\t%.1f\t%.2f\t%.2f\r\n',...
            query,Accession{i},MaxScore(i,1), TotalScore(i,1), Evalue_show(i,1), PerIdent_show(i,1));
    end
    fprintf(fid,'%s\t%s\t%.1f\t%.1f\t%.2f\t%.4f\r\n',...
        query,'HGS',max(MaxScore),sum(TotalScore),min(Evalue_show),HGS);
else %------------------------brief format --------------------
    if(isempty(header))
        return;
    end
    if(length(header)>=5)
        range=split(header{3},':');
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%.1f\t%d\t%.4f\r\n',...
            header{1},header{2},header{5},range{1},range{2},sum(MaxScore),seq_len,HGS);
    else
        fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%.1f\t%d\t%.4f\r\n',...
            header{1},'none','none','0','0',0,0,HGS);
    end
end
fclose(fid);
end
%
% %%%%%%%%%%%%%%%%%% 
%
% blastall 2.2.17   arguments:
%
%   -p  Program Name [String]
%   -d  Database [String]
%     default = nr
%   -i  Query File [File In]
%     default = stdin
%   -e  Expectation value (E) [Real]
%     default = 10.0
%   -m  alignment view options:
% 0 = pairwise,
% 1 = query-anchored showing identities,
% 2 = query-anchored no identities,
% 3 = flat query-anchored, show identities,
% 4 = flat query-anchored, no identities,
% 5 = query-anchored no identities and blunt ends,
% 6 = flat query-anchored, no identities and blunt ends,
% 7 = XML Blast output,
% 8 = tabular,
% 9 tabular with comment lines
% 10 ASN, text
% 11 ASN, binary [Integer]
%     default = 0
%     range from 0 to 11
%   -o  BLAST report Output File [File Out]  Optional
%     default = stdout
%   -F  Filter query sequence (DUST with blastn, SEG with others) [String]
%     default = T
%   -G  Cost to open a gap (-1 invokes default behavior) [Integer]
%     default = -1
%   -E  Cost to extend a gap (-1 invokes default behavior) [Integer]
%     default = -1
%   -X  X dropoff value for gapped alignment (in bits) (zero invokes default behavior)
%       blastn 30, megablast 20, tblastx 0, all others 15 [Integer]
%     default = 0
%   -I  Show GI's in deflines [T/F]
%     default = F
%   -q  Penalty for a nucleotide mismatch (blastn only) [Integer]
%     default = -3
%   -r  Reward for a nucleotide match (blastn only) [Integer]
%     default = 1
%   -v  Number of database sequences to show one-line descriptions for (V) [Integer]
%     default = 500
%   -b  Number of database sequence to show alignments for (B) [Integer]
%     default = 250
%   -f  Threshold for extending hits, default if zero
%       blastp 11, blastn 0, blastx 12, tblastn 13
%       tblastx 13, megablast 0 [Real]
%     default = 0
%   -g  Perform gapped alignment (not available with tblastx) [T/F]
%     default = T
%   -Q  Query Genetic code to use [Integer]
%     default = 1
%   -D  DB Genetic code (for tblast[nx] only) [Integer]
%     default = 1
%   -a  Number of processors to use [Integer]
%     default = 1
%   -O  SeqAlign file [File Out]  Optional
%   -J  Believe the query defline [T/F]
%     default = F
%   -M  Matrix [String]
%     default = BLOSUM62
%   -W  Word size, default if zero (blastn 11, megablast 28, all others 3) [Integer]
%     default = 0
%   -z  Effective length of the database (use zero for the real size) [Real]
%     default = 0
%   -K  Number of best hits from a region to keep (off by default, if used a value of 100 is recommended) [Integer]
%     default = 0
%   -P  0 for multiple hit, 1 for single hit (does not apply to blastn) [Integer]
%     default = 0
%   -Y  Effective length of the search space (use zero for the real size) [Real]
%     default = 0
%   -S  Query strands to search against database (for blast[nx], and tblastx)
%        3 is both, 1 is top, 2 is bottom [Integer]
%     default = 3
%   -T  Produce HTML output [T/F]
%     default = F
%   -l  Restrict search of database to list of GI's [String]  Optional
%   -U  Use lower case filtering of FASTA sequence [T/F]  Optional
%   -y  X dropoff value for ungapped extensions in bits (0.0 invokes default behavior)
%       blastn 20, megablast 10, all others 7 [Real]
%     default = 0.0
%   -Z  X dropoff value for final gapped alignment in bits (0.0 invokes default behavior)
%       blastn/megablast 50, tblastx 0, all others 25 [Integer]
%     default = 0
%   -R  PSI-TBLASTN checkpoint file [File In]  Optional
%   -n  MegaBlast search [T/F]
%     default = F
%   -L  Location on query sequence [String]  Optional
%   -A  Multiple Hits window size, default if zero (blastn/megablast 0, all others 40 [Integer]
%     default = 0
%   -w  Frame shift penalty (OOF algorithm for blastx) [Integer]
%     default = 0
%   -t  Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments. (0 invokes default behavior; a negative value disables linking.) [Integer]
%     default = 0
%   -B  Number of concatenated queries, for blastn and tblastn [Integer]  Optional
%     default = 0
%   -V  Force use of the legacy BLAST engine [T/F]  Optional
%     default = F
%   -C  Use composition-based statistics for blastp or tblastn:
%       As first character:
%       D or d: default (equivalent to T)
%       0 or F or f: no composition-based statistics
%       1 or T or t: Composition-based statistics as in NAR 29:2994-3005, 2001
%       2: Composition-based score adjustment as in Bioinformatics 21:902-911,
%           2005, conditioned on sequence properties
%       3: Composition-based score adjustment as in Bioinformatics 21:902-911,
%           2005, unconditionally
%       For programs other than tblastn, must either be absent or be D, F or 0.
%            As second character, if first character is equivalent to 1, 2, or 3:
%       U or u: unified p-value combining alignment p-value and compositional p-value in round 1 only
%  [String]
%     default = D
%   -s  Compute locally optimal Smith-Waterman alignments (This option is only
%       available for gapped tblastn.) [T/F]
%     default = F
%

%%%%%%%%%%%%%%%%%%%%%%%%
% BLASTN Program Advanced Options
%   -G  Cost to open a gap [Integer]
%     default = 5
%   -E  Cost to extend a gap [Integer]
%     default = 2
%   -q  Penalty for a mismatch in the blast portion of run [Integer]
%     default = -3
%   -r  Reward for a match in the blast portion of run [Integer]
%     default = 1
%   -e  Expectation value (E) [Real]
%     default = 10.0
%   -W  Word size, default is 11 for blastn, 3 for other programs.
%   -v  Number of one-line descriptions (V) [Integer]
%     default = 100
%   -b  Number of alignments to show (B) [Integer]
%     default = 100

%This code is KEY for calculating BLASTN.
%Weitao SUN, 2020.4.24;
function [out]=swt_fun_blastn(filefasta,db_used,fileout,format_id)
%see: http://www.molquest.com/help/2.3/programs/BlastN/parameters.html
%https://www.ncbi.nlm.nih.gov/books/NBK279684/
if(nargin<3)
    format_id=6;
end

this_cmd=strcat(['blastn -query ', filefasta,'  -evalue 8  -word_size 11  -penalty  -3  -reward  2 ',...
    '-gapopen 5  -gapextend 2 -culling_limit 1 -dust "20 64 1" -window_masker_db D:\blast\blast-2.10.0+\db\wmasker.obinary ',...
    '-db ',  db_used,' -task blastn  -outfmt ',num2str(format_id),' -out ',fileout]);
out=system(this_cmd);

%   -a  Number of processors to use [Integer],     default = 1
%see following help on blast command.
%   -r  Reward for a nucleotide match (blastn only) [Integer],      default = 1
%   -W  Word size, default if zero (blastn 11, megablast 28, all others 3) [Integer],     default = 0
%   -K  Number of best hits from a region to keep (off by default, if
%   used a value of 100 is recommended) [Integer],      default = 0
%   -q  Penalty for a nucleotide mismatch (blastn only) [Integer]
%     default = -3
%   -r  Reward for a nucleotide match (blastn only) [Integer]
%     default = 1
%   -T  Produce HTML output [T/F],    default = F
end

%weitao sun, 2004.4.23
function swt_fun_write_text(file,text,is_write_date)
if(isempty(text))
    return;
end
if(nargin<3)
    is_write_date=0;
end
date_str='';
file=swt_fun_valid_filename(file);
if(is_write_date)
    date_str=datestr(now,'yyyymmdd HH:MM:SS');
end
if(isfile(file))
    fid=fopen(file,'a');
else
    fid=fopen(file,'w');
end
if(iscell(text))
    for i=1:length(text)
        fprintf(fid,'%s\t%s\r\n',text{i},date_str);
    end
else
    fprintf(fid,'%s\t%s\r\n',text,date_str);
end
fclose(fid);
end