FasdUAS 1.101.10   ��   ��    k             l      ��  ��    ^ X Install Jon's Igor macros, and add other aliases as needed
Igor Pro signature is 	IGR0     � 	 	 �   I n s t a l l   J o n ' s   I g o r   m a c r o s ,   a n d   a d d   o t h e r   a l i a s e s   a s   n e e d e d 
 I g o r   P r o   s i g n a t u r e   i s   	 I G R 0    
  
 l     ����  r         J     
       J            m        �    F i l e   L o a d e r s   ��  m       �    H D F 5   H e l p . i h f��     ��  J           m       �    F i l e   L o a d e r s    ��   m     ! ! � " "  H D F 5 . x o p��  ��    o      ����  0 extensionslist extensionsList��  ��     # $ # l    %���� % r     & ' & b     ( ) ( o    ����  0 extensionslist extensionsList ) J     * *  + , + J     - -  . / . m     0 0 � 1 1  U t i l i t i e s /  2�� 2 m     3 3 � 4 4 ( H F S A n d P o s i x   H e l p . i h f��   ,  5�� 5 J     6 6  7 8 7 m     9 9 � : :  U t i l i t i e s 8  ;�� ; m     < < � = =  H F S A n d P o s i x . x o p��  ��   ' o      ����  0 extensionslist extensionsList��  ��   $  > ? > l     ��������  ��  ��   ?  @ A @ l     �� B C��   B : 4 Cannot run this script when Igor is already running    C � D D h   C a n n o t   r u n   t h i s   s c r i p t   w h e n   I g o r   i s   a l r e a d y   r u n n i n g A  E F E l   C G���� G Z    C H I���� H I    "�� J���� 0 appisrunning appIsRunning J  K�� K m     L L � M M  I g o r   P r o��  ��   I k   % ? N N  O P O I  % <�� Q R
�� .sysodlogaskr        TEXT Q b   % * S T S b   % ( U V U m   % & W W � X X   I g o r   i s   r u n n i n g ! V o   & '��
�� 
ret  T m   ( ) Y Y � Z Z V Q u i t   o u t   o f   I g o r   a n d   r u n   t h i s   s c r i p t   a g a i n . R �� [ \
�� 
btns [ J   + . ] ]  ^�� ^ m   + , _ _ � ` `  D o n e��   \ �� a b
�� 
dflt a m   1 2����  b �� c��
�� 
disp c m   5 6���� ��   P  d�� d L   = ?����  ��  ��  ��  ��  ��   F  e f e l  D [ g���� g I  D [�� h i
�� .sysodlogaskr        TEXT h m   D G j j � k k ` B e   P a t i e n t ,   t h i s   m a y   t a k e   u p   t o   1   m i n u t e   t o   r u n . i �� l m
�� 
btns l J   H M n n  o�� o m   H K p p � q q  C o n t i n u e��   m �� r s
�� 
dflt r m   P Q����  s �� t��
�� 
disp t m   T U���� ��  ��  ��   f  u v u l     ��������  ��  ��   v  w x w l     ��������  ��  ��   x  y z y l     �� { |��   { w q Set paths to 'Igor Pro Folder', 'Igor Extensions' folder, 'Igor Procedures' folder, and 'Igor Shared Procedures'    | � } } �   S e t   p a t h s   t o   ' I g o r   P r o   F o l d e r ' ,   ' I g o r   E x t e n s i o n s '   f o l d e r ,   ' I g o r   P r o c e d u r e s '   f o l d e r ,   a n d   ' I g o r   S h a r e d   P r o c e d u r e s ' z  ~  ~ l  \� ����� � O   \� � � � k   b� � �  � � � Q   b � � � � � r   e x � � � n   e t � � � m   p t��
�� 
ctnr � l  e p ����� � 5   e p�� ���
�� 
appf � m   i l � � � � �  I G R 0
�� kfrmID  ��  ��   � o      ���� 0 igorprofolder IgorProFolder � R      ������
�� .ascrerr ****      � ****��  ��   � k   � � � �  � � � r   � � � � � I  � ����� �
�� .sysostflalis    ��� null��   � �� ���
�� 
prmp � m   � � � � � � � 0 C h o o s e   ' I g o r   P r o   F o l d e r '��   � o      ���� 0 igorprofolder IgorProFolder �  ��� � L   � �����  ��   �  � � � l   � ��� � ���   �GA  --changed this section so xops go into Documents, not Applications
	try		set extNam to (IgorProFolder as string) & "Igor Extensions:"		set IgorExtensionsFldr to folder named extNam	on error		display dialog "Cannot find the 'Igor Extensions' Folder" buttons {"Done"} default button 1 with icon 1		return	end try
	    � � � ��     - - c h a n g e d   t h i s   s e c t i o n   s o   x o p s   g o   i n t o   D o c u m e n t s ,   n o t   A p p l i c a t i o n s 
 	 t r y  	 	 s e t   e x t N a m   t o   ( I g o r P r o F o l d e r   a s   s t r i n g )   &   " I g o r   E x t e n s i o n s : "  	 	 s e t   I g o r E x t e n s i o n s F l d r   t o   f o l d e r   n a m e d   e x t N a m  	 o n   e r r o r  	 	 d i s p l a y   d i a l o g   " C a n n o t   f i n d   t h e   ' I g o r   E x t e n s i o n s '   F o l d e r "   b u t t o n s   { " D o n e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1  	 	 r e t u r n  	 e n d   t r y  
 	 �  � � � l  � ���������  ��  ��   �  � � � l  � ��� � ���   � F @ find "Igor Shared Procedures" and save location as SharedFolder    � � � � �   f i n d   " I g o r   S h a r e d   P r o c e d u r e s "   a n d   s a v e   l o c a t i o n   a s   S h a r e d F o l d e r �  � � � Q   �� � � � � k   �B � �  � � � l  � ��� � ���   � %  first look in users local area    � � � � >   f i r s t   l o o k   i n   u s e r s   l o c a l   a r e a �  � � � r   � � � � � b   � � � � � l  � � ����� � I  � ��� � �
�� .earsffdralis        afdr � l  � � ����� � m   � ���
�� afdrdocs��  ��   � �� ���
�� 
rtyp � m   � ���
�� 
TEXT��  ��  ��   � m   � � � � � � � � W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s : I g o r   S h a r e d   P r o c e d u r e s : � o      ���� 0 
sharedname 
SharedName �  � � � r   � � � � � 5   � ��� ���
�� 
cfol � o   � ����� 0 
sharedname 
SharedName
�� kfrmname � o      ���� 0 sharedfolder SharedFolder �  � � � r   � � � � � m   � � � � � � � @ i n s t a l l i n g   i n t o   D o c u m e n t s   F o l d e r � o      ���� 0 outcome   �  � � � r   � � � � � m   � ���
�� boovtrue � o      ���� 0 
usinglocal 
usingLocal �  � � � l  � ���������  ��  ��   �  � � � r   � � � � � 5   � ��� ���
�� 
cfol � l  � � ����� � b   � � � � � l  � � ����� � I  � ��� � �
�� .earsffdralis        afdr � l  � � ����� � m   � ���
�� afdrdocs��  ��   � �� ���
�� 
rtyp � m   � ���
�� 
TEXT��  ��  ��   � m   � � � � � � � d W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s :��  ��  
�� kfrmname � o      ����  0 userprocedures UserProcedures �  ��� � Z   �B � ����� � H   � � � � l  � � ����� � I  � ��� ���
�� .coredoexbool        obj  � n   � � � � � 4   � ��� �
�� 
cfol � m   � � � � � � �  L o c a l   P a c k a g e s � o   � �����  0 userprocedures UserProcedures��  ��  ��   � k   �> � �  � � � I  � ���� �
�� .corecrel****      � null��   � �� � �
�� 
kocl � m   ���
�� 
cfol � �� � �
�� 
insh � o  ����  0 userprocedures UserProcedures � �� ��
�� 
prdt � K  
 � � �~ � �
�~ 
pnam � m   � � � � �  L o c a l   P a c k a g e s � �} ��|
�} 
labi � m  �{�{ �|  �   �  � � � r  !5 � � � m  !$�z�z  � n       � � � 1  04�y
�y 
labi � n  $0 �  � 4  )0�x
�x 
cfol m  ,/ � � U s e r s : t i s c h l e r : D o c u m e n t s : W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s : L o c a l   P a c k a g e s  1  $)�w
�w 
sdsk � �v n 6> I  7>�u�t�u 0 makeaboutfile makeAboutFile �s o  7:�r�r  0 userprocedures UserProcedures�s  �t    f  67�v  ��  ��  ��   � R      �q�p�o
�q .ascrerr ****      � ****�p  �o   � k  J�		 

 l JJ�n�n   7 1 not in users local area, look in Igor Pro Folder    � b   n o t   i n   u s e r s   l o c a l   a r e a ,   l o o k   i n   I g o r   P r o   F o l d e r  Q  J� k  Mz  r  M\ b  MX l MT�m�l c  MT o  MP�k�k 0 igorprofolder IgorProFolder m  PS�j
�j 
TEXT�m  �l   m  TW � N U s e r   P r o c e d u r e s : I g o r   S h a r e d   P r o c e d u r e s : o      �i�i 0 
sharedname 
SharedName  !  r  ]l"#" 5  ]h�h$�g
�h 
cfol$ o  ad�f�f 0 
sharedname 
SharedName
�g kfrmname# o      �e�e 0 sharedfolder SharedFolder! %&% r  mt'(' m  mp)) �** F i n s t a l l i n g   i n t o   A p p l i c a t i o n s   F o l d e r( o      �d�d 0 outcome  & +�c+ r  uz,-, m  uv�b
�b boovfals- o      �a�a 0 
usinglocal 
usingLocal�c   R      �`�_�^
�` .ascrerr ****      � ****�_  �^   k  ��.. /0/ l ���]12�]  1 #  cannot find it, ask the user   2 �33 :   c a n n o t   f i n d   i t ,   a s k   t h e   u s e r0 454 r  ��676 I ���\�[8
�\ .sysostflalis    ��� null�[  8 �Z9�Y
�Z 
prmp9 m  ��:: �;; > C h o o s e   ' I g o r   S h a r e d   P r o c e d u r e s '�Y  7 o      �X�X 0 sharedfolder SharedFolder5 <�W< l ���V=>�V  =   set SharedFolder to none   > �?? 2   s e t   S h a r e d F o l d e r   t o   n o n e�W   @A@ l ���UBC�U  B U Oset SharedFolder to choose folder with prompt "Choose 'Igor Shared Procedures'"   C �DD � s e t   S h a r e d F o l d e r   t o   c h o o s e   f o l d e r   w i t h   p r o m p t   " C h o o s e   ' I g o r   S h a r e d   P r o c e d u r e s ' "A E�TE l ���SFG�S  F   set SharedFolder to none   G �HH 2   s e t   S h a r e d F o l d e r   t o   n o n e�T   � IJI l ���RKL�R  K n h IgorProcedures points to the folder that will hold 'always alias' (either in Documents or Applications)   L �MM �   I g o r P r o c e d u r e s   p o i n t s   t o   t h e   f o l d e r   t h a t   w i l l   h o l d   ' a l w a y s   a l i a s '   ( e i t h e r   i n   D o c u m e n t s   o r   A p p l i c a t i o n s )J NON r  ��PQP n  ��RSR 1  ���Q
�Q 
pareS l ��T�P�OT n  ��UVU 1  ���N
�N 
pareV o  ���M�M 0 sharedfolder SharedFolder�P  �O  Q o      �L�L 0 	topfolder 	topFolderO WXW r  ��YZY 4  ���K[
�K 
alis[ l ��\�J�I\ b  ��]^] l ��_�H�G_ c  ��`a` o  ���F�F 0 	topfolder 	topFoldera m  ���E
�E 
TEXT�H  �G  ^ m  ��bb �cc  I g o r   P r o c e d u r e s�J  �I  Z o      �D�D  0 igorprocedures IgorProceduresX ded l ���C�B�A�C  �B  �A  e fgf l ���@hi�@  h p j IgorExtensionsFldr points to the folder that will hold xop aliases, (either in Documents or Applications)   i �jj �   I g o r E x t e n s i o n s F l d r   p o i n t s   t o   t h e   f o l d e r   t h a t   w i l l   h o l d   x o p   a l i a s e s ,   ( e i t h e r   i n   D o c u m e n t s   o r   A p p l i c a t i o n s )g k�?k r  ��lml 5  ���>n�=
�> 
cfoln l ��o�<�;o b  ��pqp l ��r�:�9r c  ��sts o  ���8�8 0 	topfolder 	topFoldert m  ���7
�7 
TEXT�:  �9  q m  ��uu �vv  I g o r   E x t e n s i o n s�<  �;  
�= kfrmnamem o      �6�6 (0 igorextensionsfldr IgorExtensionsFldr�?   � m   \ _ww�                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  ��  ��    xyx l     �5�4�3�5  �4  �3  y z{z l     �2�1�0�2  �1  �0  { |}| l     �/~�/  ~ ) # Set path to 'always' source folder    ��� F   S e t   p a t h   t o   ' a l w a y s '   s o u r c e   f o l d e r} ��� l ���.�-� Q  ����� r  ����� 4  ���,�
�, 
alis� l ����+�*� b  ����� l ����)�(� c  ����� o  ���'�' 0 sharedfolder SharedFolder� m  ���&
�& 
TEXT�)  �(  � m  ���� ���  a l w a y s�+  �*  � o      �%�% $0 alwaysfolderorig alwaysFolderOrig� R      �$�#�"
�$ .ascrerr ****      � ****�#  �"  � k  ��� ��� I �
�!��
�! .sysodlogaskr        TEXT� m  ���� ��� f C a n n o t   f i n d   f o l d e r   ' : I g o r   S h a r e d   P r o c e d u r e s : a l w a y s '� � ��
�  
btns� J  ���� ��� m  ���� ���  D o n e�  � ���
� 
dflt� m  � �� � ���
� 
disp� m  �� �  � ��� L  ��  �  �.  �-  � ��� l     ����  �  �  � ��� l     ����  �  �  � ��� l +���� r  +��� b  '��� b  ��� b  ��� b  ��� b  ��� b  ��� o  �� 0 outcome  � o  �
� 
ret � m  �� ���  - - - - - - - - - - - - -� o  �
� 
ret � m  �� ��� ( u s i n g   I g o r P r o F o l d e r :� o  �
� 
ret � l &���
� c  &��� o  "�	�	 0 igorprofolder IgorProFolder� m  "%�
� 
TEXT�  �
  � o      �� 0 outcome  �  �  � ��� l     ����  � P Jset outcome to "using IgorProFolder:" & return & (IgorProFolder as string)   � ��� � s e t   o u t c o m e   t o   " u s i n g   I g o r P r o F o l d e r : "   &   r e t u r n   &   ( I g o r P r o F o l d e r   a s   s t r i n g )� ��� l     ����  � � �display dialog "using IgorProFolder:" & return & return & (IgorProFolder as string) buttons {"Continue"} default button 1 with icon 1   � ���
 d i s p l a y   d i a l o g   " u s i n g   I g o r P r o F o l d e r : "   &   r e t u r n   &   r e t u r n   &   ( I g o r P r o F o l d e r   a s   s t r i n g )   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1� ��� l     ����  � } wdisplay dialog "IgorProcedures" & return & (IgorProcedures as string) buttons {"Continue"} default button 1 with icon 1   � ��� � d i s p l a y   d i a l o g   " I g o r P r o c e d u r e s "   &   r e t u r n   &   ( I g o r P r o c e d u r e s   a s   s t r i n g )   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1� ��� l     ����  � � {display dialog "alwaysFolderOrig" & return & (alwaysFolderOrig as string) buttons {"Continue"} default button 1 with icon 1   � ��� � d i s p l a y   d i a l o g   " a l w a y s F o l d e r O r i g "   &   r e t u r n   &   ( a l w a y s F o l d e r O r i g   a s   s t r i n g )   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1� ��� l     ����  � � display dialog "IgorExtensionsFldr" & return & (IgorExtensionsFldr as string) buttons {"Continue"} default button 1 with icon 1   � ��� � d i s p l a y   d i a l o g   " I g o r E x t e n s i o n s F l d r "   &   r e t u r n   &   ( I g o r E x t e n s i o n s F l d r   a s   s t r i n g )   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1� ��� l     ����  � d ^display dialog "usingLocal is " & usingLocal buttons {"Continue"} default button 1 with icon 1   � ��� � d i s p l a y   d i a l o g   " u s i n g L o c a l   i s   "   &   u s i n g L o c a l   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1� ��� l     � �����   ��  ��  � ��� l     ��������  ��  ��  � ��� l     ������  � W Q Check that the alias to 'always' in 'Igor Procedures' is present, if not make it   � ��� �   C h e c k   t h a t   t h e   a l i a s   t o   ' a l w a y s '   i n   ' I g o r   P r o c e d u r e s '   i s   p r e s e n t ,   i f   n o t   m a k e   i t� ��� l ,������� Z  ,������� H  ,8�� l ,7������ I  ,7�������  0 isaliascorrect isAliasCorrect� ��� o  -0����  0 igorprocedures IgorProcedures� ���� o  03���� $0 alwaysfolderorig alwaysFolderOrig��  ��  ��  ��  � k  ;��� ��� O  ;���� Q  A����� r  Da��� I D]�����
�� .corecrel****      � null��  � ����
�� 
kocl� m  HK��
�� 
alia� ��� 
�� 
insh� o  NQ����  0 igorprocedures IgorProcedures  ����
�� 
to   o  TW���� $0 alwaysfolderorig alwaysFolderOrig��  � o      ���� 0 alwaysalias alwaysAlias� R      ������
�� .ascrerr ****      � ****��  ��  � k  i�  I i���
�� .sysodlogaskr        TEXT m  il � N C o u l d   n o t   m a k e   a l i a s   t o   ' a l w a y s '   f o l d e r ��	

�� 
btns	 J  mr �� m  mp �  D o n e��  
 ��
�� 
dflt m  uv����  ����
�� 
disp m  yz���� ��   �� L  ������  ��  � m  ;>�                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  �  r  �� b  �� b  �� b  �� b  �� o  ������ 0 outcome   o  ����
�� 
ret  m  ��   �!!  - - - - - - - - - - - - - o  ����
�� 
ret  m  ��"" �## : m a d e   a l i a s   t o   ' a l w a y s '   f o l d e r o      ���� 0 outcome   $��$ l ����%&��  % f `display dialog "made alias to 'always' folder" buttons {"Continue"} default button 1 with icon 1   & �'' � d i s p l a y   d i a l o g   " m a d e   a l i a s   t o   ' a l w a y s '   f o l d e r "   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1��  ��  � k  ��(( )*) r  ��+,+ b  ��-.- b  ��/0/ b  ��121 b  ��343 b  ��565 b  ��787 o  ������ 0 outcome  8 o  ����
�� 
ret 6 m  ��99 �::  - - - - - - - - - - - - -4 o  ����
�� 
ret 2 m  ��;; �<< * ' a l w a y s   a l i a s '   i s   O K ,0 o  ����
�� 
ret . m  ��== �>>   n o t h i n g   c h a n g e d ., o      ���� 0 outcome  * ?��? l ����@A��  @ | vdisplay dialog "'always alias' is OK," & return & "nothing changed." buttons {"Continue"} default button 1 with icon 1   A �BB � d i s p l a y   d i a l o g   " ' a l w a y s   a l i a s '   i s   O K , "   &   r e t u r n   &   " n o t h i n g   c h a n g e d . "   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1��  ��  ��  � CDC l     ��������  ��  ��  D EFE l     ��������  ��  ��  F GHG l     ��������  ��  ��  H IJI l     ��KL��  K N H Check that the aliases to all of the xops are present, if not make them   L �MM �   C h e c k   t h a t   t h e   a l i a s e s   t o   a l l   o f   t h e   x o p s   a r e   p r e s e n t ,   i f   n o t   m a k e   t h e mJ NON l ��P����P r  ��QRQ m  ��SS �TT  R o      ���� $0 unchangedaliases unchangedAliases��  ��  O UVU l ��W����W r  ��XYX m  ��ZZ �[[  Y o      ����  0 changedaliases changedAliases��  ��  V \]\ l �^����^ X  �_��`_ k  ��aa bcb r  ��ded n  ��fgf 4  ����h
�� 
cobjh m  ������ g o  ������ 0 
iteminlist 
itemInListe o      ���� 0 thepath thePathc iji r  ��klk n  ��mnm 4  ����o
�� 
cobjo m  ������ n o  ������ 0 
iteminlist 
itemInListl o      ���� 0 thename theNamej pqp l ����������  ��  ��  q rsr Q  �xtuvt k  �(ww xyx r  �z{z b  �|}| b  �	~~ b  ���� b  ���� l �������� c  ����� o  ������ 0 igorprofolder IgorProFolder� m  ����
�� 
TEXT��  ��  � m  � �� ���   M o r e   E x t e n s i o n s :� o  ���� 0 thepath thePath m  �� ���  :} o  	���� 0 thename theName{ o      ���� 0 origname origNamey ���� O  (��� r  '��� 5  #�����
�� 
file� o  ���� 0 origname origName
�� kfrmname� o      ���� 0 original  � m  ���                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  ��  u R      ������
�� .ascrerr ****      � ****��  ��  v k  0x�� ��� I 0O����
�� .sysodlogaskr        TEXT� b  0;��� b  07��� m  03�� ��� H C a n n o t   f i n d   t h e   o r i g n a l   v e r s i o n   o f   '� o  36���� 0 thename theName� m  7:�� ���  '� ����
�� 
btns� J  <A�� ���� m  <?�� ���  C o n t i n u e��  � ����
�� 
dflt� m  DE���� � �����
�� 
disp� m  HI���� ��  � ��� I Pu����
�� .sysodlogaskr        TEXT� b  Pa��� b  P]��� b  PY��� b  PU��� m  PS�� ���  w a s   l o o k i n g   i n� o  ST��
�� 
ret � m  UX�� ���  '� o  Y\���� 0 origname origName� m  ]`�� ���  '� ����
�� 
btns� J  bg�� ���� m  be�� ���  D o n e��  � ����
�� 
dflt� m  jk���� � �����
�� 
disp� m  no���� ��  � ���� L  vx����  ��  s ��� l yy��������  ��  ��  � ���� Z  y������� H  y��� l y������� I  y��������  0 isaliascorrect isAliasCorrect� ��� o  z}���� (0 igorextensionsfldr IgorExtensionsFldr� ��� o  }��~�~ 0 original  �  ��  ��  ��  � O  ����� k  ���� ��� r  ����� b  ����� b  ����� b  ����� b  ����� o  ���}�}  0 changedaliases changedAliases� o  ���|
�| 
ret � o  ���{�{ 0 thepath thePath� m  ���� ���  :� o  ���z�z 0 thename theName� o      �y�y  0 changedaliases changedAliases� ��x� Q  ������ r  ����� I ���w�v�
�w .corecrel****      � null�v  � �u��
�u 
kocl� m  ���t
�t 
alia� �s��
�s 
insh� o  ���r�r (0 igorextensionsfldr IgorExtensionsFldr� �q��p
�q 
to  � o  ���o�o 0 original  �p  � o      �n�n 0 repeatalias repeatAlias� R      �m�l�k
�m .ascrerr ****      � ****�l  �k  � k  ���� ��� I ���j��
�j .sysodlogaskr        TEXT� b  ����� b  ����� m  ���� ��� 2 C o u l d   n o t   m a k e   a l i a s   t o   '� o  ���i�i 0 thename theName� m  ���� ���  '� �h��
�h 
btns� J  ���� ��g� m  ���� ���  D o n e�g  � �f 
�f 
dflt  m  ���e�e  �d�c
�d 
disp m  ���b�b �c  � �a L  ���`�`  �a  �x  � m  ���                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  ��  � r  �� b  �� b  ��	
	 o  ���_�_ $0 unchangedaliases unchangedAliases
 o  ���^
�^ 
ret  o  ���]�] 0 thename theName o      �\�\ $0 unchangedaliases unchangedAliases��  �� 0 
iteminlist 
itemInList` o  ���[�[  0 extensionslist extensionsList��  ��  ]  l     �Z�Y�X�Z  �Y  �X    l ,�W�V Z  ,�U�T ?  n   1  �S
�S 
leng o  �R�R $0 unchangedaliases unchangedAliases m  �Q�Q   r  ( b  $ b    b   b   b   !  o  �P�P 0 outcome  ! o  �O
�O 
ret  m  "" �##  - - - - - - - - - - - - - o  �N
�N 
ret  m  $$ �%% 4 D i d   n o t   c h a n g e   a l i a s e s   t o : o   #�M�M $0 unchangedaliases unchangedAliases o      �L�L 0 outcome  �U  �T  �W  �V   &'& l -T(�K�J( Z  -T)*�I�H) ? -6+,+ n  -4-.- 1  04�G
�G 
leng. o  -0�F�F  0 changedaliases changedAliases, m  45�E�E  * r  9P/0/ b  9L121 b  9H343 b  9D565 b  9B787 b  9>9:9 o  9<�D�D 0 outcome  : o  <=�C
�C 
ret 8 m  >A;; �<<  - - - - - - - - - - - - -6 o  BC�B
�B 
ret 4 m  DG== �>> " C r e a t e d   a l i a s   t o :2 o  HK�A�A  0 changedaliases changedAliases0 o      �@�@ 0 outcome  �I  �H  �K  �J  ' ?@? l     �?�>�=�?  �>  �=  @ ABA l     �<�;�:�<  �;  �:  B CDC l USE�9�8E Z  USFG�7�6F = UZHIH o  UX�5�5 0 
usinglocal 
usingLocalI m  XY�4
�4 boovtrueG k  ]OJJ KLK O ]{MNM r  czOPO 5  cv�3Q�2
�3 
cfolQ l grR�1�0R b  grSTS l gnU�/�.U c  gnVWV o  gj�-�- 0 igorprofolder IgorProFolderW m  jm�,
�, 
TEXT�/  �.  T m  nqXX �YY  I g o r   E x t e n s i o n s�1  �0  
�2 kfrmnameP o      �+�+ .0 igorextensionsfldrapp IgorExtensionsFldrAppN m  ]`ZZ�                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  L [\[ r  |�]^] m  |__ �``  ^ o      �*�* 0 duplicatexops duplicateXOPs\ aba X  �c�)dc k  � ee fgf r  ��hih n  ��jkj 4  ���(l
�( 
cobjl m  ���'�' k o  ���&�& 0 
iteminlist 
itemInListi o      �%�% 0 thepath thePathg mnm r  ��opo n  ��qrq 4  ���$s
�$ 
cobjs m  ���#�# r o  ���"�" 0 
iteminlist 
itemInListp o      �!�! 0 thename theNamen tut r  ��vwv b  ��xyx b  ��z{z b  ��|}| b  ��~~ l ���� �� c  ����� o  ���� 0 igorprofolder IgorProFolder� m  ���
� 
TEXT�   �   m  ���� ���   M o r e   E x t e n s i o n s :} o  ���� 0 thepath thePath{ m  ���� ���  :y o  ���� 0 thename theNamew o      �� 0 origname origNameu ��� O ����� r  ����� 5  �����
� 
file� o  ���� 0 origname origName
� kfrmname� o      �� 0 original  � m  �����                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  � ��� Z  � ����� I  ������  0 isaliascorrect isAliasCorrect� ��� o  ���� .0 igorextensionsfldrapp IgorExtensionsFldrApp� ��� o  ���� 0 original  �  �  � r  ����� b  ����� b  ����� o  ���� 0 duplicatexops duplicateXOPs� o  ���
� 
ret � o  ���� 0 thename theName� o      �
�
 0 duplicatexops duplicateXOPs�  �  �  �) 0 
iteminlist 
itemInListd o  ���	�	  0 extensionslist extensionsListb ��� Z  O����� ? ��� n  ��� 1  	�
� 
leng� o  	�� 0 duplicatexops duplicateXOPs� m  ��  � k  K�� ��� I ��� 
� .sysobeepnull��� ��� long�  �   � ���� I K����
�� .sysodlogaskr        TEXT� b  7��� b  3��� b  1��� b  /��� b  '��� b  #��� b  !��� b  ��� m  �� ��� * F o u n d   D u p l i c a t e   X O P s :� o  ���� 0 duplicatexops duplicateXOPs� o   ��
�� 
ret � o  !"��
�� 
ret � m  #&�� ���  i n  � l '.������ c  '.��� o  '*���� .0 igorextensionsfldrapp IgorExtensionsFldrApp� m  *-��
�� 
TEXT��  ��  � o  /0��
�� 
ret � o  12��
�� 
ret � m  36�� ��� @ Y o u   s h o u l d   p r o b a b l y   d e l e t e   t h i s !� ����
�� 
btns� J  8=�� ���� m  8;�� ���  C o n t i n u e��  � ����
�� 
dflt� m  @A���� � �����
�� 
disp� m  DE���� ��  ��  �  �  �  �7  �6  �9  �8  D ��� l     ��������  ��  ��  � ��� l     ��������  ��  ��  � ��� l     ��������  ��  ��  � ��� l     ������  � M Gset outcome to outcome & return & "-------------" & return & "Finished"   � ��� � s e t   o u t c o m e   t o   o u t c o m e   &   r e t u r n   &   " - - - - - - - - - - - - - "   &   r e t u r n   &   " F i n i s h e d "� ��� l Tk������ I Tk����
�� .sysodlogaskr        TEXT� o  TW���� 0 outcome  � ����
�� 
btns� J  X]�� ���� m  X[�� ���  D o n e��  � ����
�� 
dflt� m  `a���� � �����
�� 
disp� m  de���� ��  ��  ��  � ��� l ln������ L  ln����  ��  ��  � ��� l     ������  �   DONE   � ��� 
   D O N E� ��� l     ��������  ��  ��  � ��� l     ��������  ��  ��  � ��� i     ��� I      �������  0 isaliascorrect isAliasCorrect� ��� o      ���� 0 fldr  � ���� o      ���� 0 original  ��  ��  � k     G�� ��� l     ������  � Q K in the folder fldr, check to see if there is already an alias to original	   � ��� �   i n   t h e   f o l d e r   f l d r ,   c h e c k   t o   s e e   i f   t h e r e   i s   a l r e a d y   a n   a l i a s   t o   o r i g i n a l 	� ��� r        n     2   ��
�� 
file o     ���� 0 fldr   o      ���� 0 filelist fileList�  O    D X   
 C��	 Z    >
����
 =   ! c     l   ���� n     1    ��
�� 
kind o    ���� 0 
iteminlist 
itemInList��  ��   m    ��
�� 
TEXT m      � 
 A l i a s k   $ :  r   $ ) n   $ ' 1   % '��
�� 
orig o   $ %���� 0 
iteminlist 
itemInList o      ���� 0 
sourcefile 
sourceFile �� Z   * :���� =  * 1  l  * -!����! c   * -"#" o   * +���� 0 
sourcefile 
sourceFile# m   + ,��
�� 
TEXT��  ��    l  - 0$����$ c   - 0%&% o   - .���� 0 original  & m   . /��
�� 
TEXT��  ��   L   4 6'' m   4 5��
�� boovtrue��  ��  ��  ��  ��  �� 0 
iteminlist 
itemInList	 o    ���� 0 filelist fileList m    ((�                                                                                  MACS  alis    r  Macintosh HD               Ũ�H+  ?�8
Finder.app                                                     @"�Ƙ�        ����  	                CoreServices    ũz      ƘK�    ?�8?��?��  3Macintosh HD:System:Library:CoreServices:Finder.app    
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��   )��) L   E G** m   E F��
�� boovfals��  � +,+ l     ��������  ��  ��  , -.- l     ��������  ��  ��  . /0/ i    121 I      ��3���� 0 appisrunning appIsRunning3 4��4 o      ���� 0 appname AppName��  ��  2 O    565 E    787 l   	9����9 n    	:;: 1    	��
�� 
pnam; 2   ��
�� 
prcs��  ��  8 o   	 
���� 0 appname AppName6 m     <<�                                                                                  sevs  alis    �  Macintosh HD               Ũ�H+  ?�8System Events.app                                              @���7��        ����  	                CoreServices    ũz      �8'7    ?�8?��?��  :Macintosh HD:System:Library:CoreServices:System Events.app  $  S y s t e m   E v e n t s . a p p    M a c i n t o s h   H D  -System/Library/CoreServices/System Events.app   / ��  0 =>= l     ��������  ��  ��  > ?@? l     ��������  ��  ��  @ ABA i    CDC I      ��E���� 0 makeaboutfile makeAboutFileE F��F o      ����  0 userprocedures UserProcedures��  ��  D O     |GHG k    {II JKJ I   	������
�� .miscactvnull��� ��� null��  ��  K LML I  
 ����N
�� .corecrel****      � null��  N ��O��
�� 
koclO m    ��
�� 
docu��  M PQP r    RSR m    TT �UU� P u t   a n y   I g o r   p a c k a g e s   i n   t h e   t o p   l e v e l   o f   t h i s   f o l d e r   t h a t   y o u   w a n t   a v a i l a b l e   f o r   l o a d i n g .     T h e y   w i l l   a u t o m a t i c a l l y   b e   a v a i l a b l e   f o r   i n c l u d i n g   b y   g o i n g   t o   t h e   I g o r   M e n u       " A n a l y s i s : P a c k a g e s : A d d   L o c a l   U s e r   P a c k a g e s . . . "S o      ���� 0 body  Q VWV r    XYX b    Z[Z b    \]\ b    ^_^ o    ���� 0 body  _ o    ��
�� 
ret ] o    ��
�� 
ret [ m    `` �aa z I f   i n   t h e   f i r s t   5 0   l i n e s   o f   t h e   f i l e   t h e r e   i s   a   l i n e   s u c h   a s :Y o      ���� 0 body  W bcb r     +ded b     )fgf b     'hih b     %jkj b     #lml o     !���� 0 body  m o   ! "��
�� 
ret k o   # $��
�� 
ret i m   % &nn �oo \ # r e q u i r e d P a c k a g e s   " H D F 5 i m a g e s ; m i c r o G e o m e t r y N ; "g o   ' (��
�� 
ret e o      ���� 0 body  c pqp r   , 3rsr b   , 1tut b   , /vwv o   , -���� 0 body  w o   - .�
� 
ret u m   / 0xx �yy� t h e n   t h i s   i p f   f i l e   w i l l   o n l y   b e   l o a d e d   i f   a l l   o f   t h e   i p f   f i l e s   i n   t h e   l i s t   a r e   a l r e a d y   l o a d e d .     I n   t h i s   c a s e ,   t h e   i p f   f i l e   w i l l   o n l y   l o a d   i f   b o t h   H D F 5 i m a g e s   a n d   m i c r o G e o m e t r y N   h a v e   a l r e a d y   b e e n   l o a d e d .s o      �~�~ 0 body  q z{z r   4 ;|}| b   4 9~~ b   4 7��� o   4 5�}�} 0 body  � o   5 6�|
�| 
ret  m   7 8�� ���� P u t   a n y   I g o r   p a c k a g e s   i n   t h e   t o p   l e v e l   o f   t h i s   f o l d e r   t h a t   y o u   w a n t   a v a i l a b l e   f o r   l o a d i n g .     T h e y   w i l l   a u t o m a t i c a l l y   b e   a v a i l a b l e   f o r   i n c l u d i n g   b y   g o i n g   t o   t h e   I g o r   M e n u       " A n a l y s i s : P a c k a g e s : A d d   L o c a l   U s e r   P a c k a g e s . . . "} o      �{�{ 0 body  { ��� r   < C��� b   < A��� b   < ?��� o   < =�z�z 0 body  � o   = >�y
�y 
ret � m   ? @�� ��� z I f   i n   t h e   f i r s t   5 0   l i n e s   o f   t h e   f i l e   t h e r e   i s   a   l i n e   s u c h   a s :� o      �x�x 0 body  � ��� r   D K��� b   D I��� b   D G��� o   D E�w�w 0 body  � o   E F�v
�v 
ret � m   G H�� ��� \ # r e q u i r e d P a c k a g e s   " H D F 5 i m a g e s ; m i c r o G e o m e t r y N ; "� o      �u�u 0 body  � ��� r   L S��� b   L Q��� b   L O��� o   L M�t�t 0 body  � o   M N�s
�s 
ret � m   O P�� ���� t h e n   t h i s   i p f   f i l e   w i l l   o n l y   b e   l o a d e d   i f   a l l   o f   t h e   i p f   f i l e s   i n   t h e   l i s t   a r e   a l r e a d y   l o a d e d .     I n   t h i s   c a s e ,   t h e   i p f   f i l e   w i l l   o n l y   l o a d   i f   b o t h   H D F 5 i m a g e s   a n d   m i c r o G e o m e t r y N   h a v e   a l r e a d y   b e e n   l o a d e d .� o      �r�r 0 body  � ��� r   T ^��� c   T W��� o   T U�q�q 0 body  � m   U V�p
�p 
ctxt� n      ��� m   [ ]�o
�o 
ctxt� 4   W [�n�
�n 
docu� m   Y Z�m�m � ��� r   _ h��� b   _ f��� l  _ b��l�k� c   _ b��� o   _ `�j�j  0 userprocedures UserProcedures� m   ` a�i
�i 
TEXT�l  �k  � m   b e�� ��� T L o c a l   P a c k a g e s : A b o u t   L o c a l   P a c k a g e s   F o l d e r� o      �h�h 0 	aboutfile 	aboutFile� ��� I  i u�g��
�g .coresavenull���    obj � 4   i m�f�
�f 
docu� m   k l�e�e � �d��c
�d 
kfil� o   p q�b�b 0 	aboutfile 	aboutFile�c  � ��a� I  v {�`�_�^
�` .aevtquitnull��� ��� null�_  �^  �a  H m     ���                                                                                  ttxt  alis    T  Macintosh HD               Ũ�H+  ?�@TextEdit.app                                                   @�p�k!        ����  	                Applications    ũz      �kgd    ?�@  &Macintosh HD:Applications:TextEdit.app    T e x t E d i t . a p p    M a c i n t o s h   H D  Applications/TextEdit.app   / ��  B ��� l     �]�\�[�]  �\  �[  � ��� l     �Z���Z  �  tell application "Finder"   � ��� 2 t e l l   a p p l i c a t i o n   " F i n d e r "� ��� l     �Y���Y  � G A    set theFile to "/Users/edstockly/Documents/2-CoverHighlights"   � ��� �         s e t   t h e F i l e   t o   " / U s e r s / e d s t o c k l y / D o c u m e n t s / 2 - C o v e r H i g h l i g h t s "� ��� l     �X���X  � . (    set theFile to theFile as POSIX file   � ��� P         s e t   t h e F i l e   t o   t h e F i l e   a s   P O S I X   f i l e� ��� l     �W���W  � D >    set aliasFile to make new alias file at desktop to theFile   � ��� |         s e t   a l i a s F i l e   t o   m a k e   n e w   a l i a s   f i l e   a t   d e s k t o p   t o   t h e F i l e� ��� l     �V���V  � 7 1    set the name of aliasFile to "New Alias File"   � ��� b         s e t   t h e   n a m e   o f   a l i a s F i l e   t o   " N e w   A l i a s   F i l e "� ��� l     �U���U  �  end tell   � ���  e n d   t e l l� ��� l     �T�S�R�T  �S  �R  � ��Q� l     �P�O�N�P  �O  �N  �Q       "�M�����������L������Z 9 <���_�K�J�I�H�G�F�E�D�C�M  �  �B�A�@�?�>�=�<�;�:�9�8�7�6�5�4�3�2�1�0�/�.�-�,�+�*�)�(�'�&�%�$�#�B  0 isaliascorrect isAliasCorrect�A 0 appisrunning appIsRunning�@ 0 makeaboutfile makeAboutFile
�? .aevtoappnull  �   � ****�>  0 extensionslist extensionsList�= 0 igorprofolder IgorProFolder�< 0 
sharedname 
SharedName�; 0 sharedfolder SharedFolder�: 0 outcome  �9 0 
usinglocal 
usingLocal�8  0 userprocedures UserProcedures�7 0 	topfolder 	topFolder�6  0 igorprocedures IgorProcedures�5 (0 igorextensionsfldr IgorExtensionsFldr�4 $0 alwaysfolderorig alwaysFolderOrig�3 $0 unchangedaliases unchangedAliases�2  0 changedaliases changedAliases�1 0 thepath thePath�0 0 thename theName�/ 0 origname origName�. 0 original  �- .0 igorextensionsfldrapp IgorExtensionsFldrApp�, 0 duplicatexops duplicateXOPs�+  �*  �)  �(  �'  �&  �%  �$  �#  � �"��!� ����"  0 isaliascorrect isAliasCorrect�! ��� �  ��� 0 fldr  � 0 original  �   � ������ 0 fldr  � 0 original  � 0 filelist fileList� 0 
iteminlist 
itemInList� 0 
sourcefile 
sourceFile� 	�(������
� 
file
� 
kocl
� 
cobj
� .corecnte****       ****
� 
kind
� 
TEXT
� 
orig� H��-E�O� ; 8�[��l kh ��,�&�  ��,E�O��&��&  eY hY h[OY��UOf� �2������ 0 appisrunning appIsRunning� ��� �  �
�
 0 appname AppName�  � �	�	 0 appname AppName� <��
� 
prcs
� 
pnam� � 	*�-�,�U� �D������ 0 makeaboutfile makeAboutFile� ��� �  ��  0 userprocedures UserProcedures�  � � �����   0 userprocedures UserProcedures�� 0 body  �� 0 	aboutfile 	aboutFile� ���������T��`nx���������������
�� .miscactvnull��� ��� null
�� 
kocl
�� 
docu
�� .corecrel****      � null
�� 
ret 
�� 
ctxt
�� 
TEXT
�� 
kfil
�� .coresavenull���    obj 
�� .aevtquitnull��� ��� null� }� y*j O*��l O�E�O��%�%�%E�O��%�%�%�%E�O��%�%E�O��%�%E�O��%�%E�O��%�%E�O��%�%E�O��&*�k/�-FO��&a %E�O*�k/a �l O*j U� �����������
�� .aevtoappnull  �   � ****� k    n    
  #  E  e  ~ � � � N		 U

 \  & C � �����  ��  ��  � ���� 0 
iteminlist 
itemInList� �    !�� 0 3 9 < L�� W�� Y�� _�������� j pw�� ������������� ����������� ��������� ����� ��� ����������� �����������):������b��u����������������� "9;=S��Z���������������������������������"$;=X��_�������������  0 extensionslist extensionsList�� 0 appisrunning appIsRunning
�� 
ret 
�� 
btns
�� 
dflt
�� 
disp�� 
�� .sysodlogaskr        TEXT
�� 
appf
�� kfrmID  
�� 
ctnr�� 0 igorprofolder IgorProFolder��  ��  
�� 
prmp
�� .sysostflalis    ��� null
�� afdrdocs
�� 
rtyp
�� 
TEXT
�� .earsffdralis        afdr�� 0 
sharedname 
SharedName
�� 
cfol
�� kfrmname�� 0 sharedfolder SharedFolder�� 0 outcome  �� 0 
usinglocal 
usingLocal��  0 userprocedures UserProcedures
�� .coredoexbool        obj 
�� 
kocl
�� 
insh
�� 
prdt
�� 
pnam
�� 
labi�� 
�� .corecrel****      � null
�� 
sdsk�� 0 makeaboutfile makeAboutFile
�� 
pare�� 0 	topfolder 	topFolder
�� 
alis��  0 igorprocedures IgorProcedures�� (0 igorextensionsfldr IgorExtensionsFldr�� $0 alwaysfolderorig alwaysFolderOrig��  0 isaliascorrect isAliasCorrect
�� 
alia
�� 
to  �� 0 alwaysalias alwaysAlias�� $0 unchangedaliases unchangedAliases��  0 changedaliases changedAliases
�� 
cobj
�� .corecnte****       ****�� 0 thepath thePath�� 0 thename theName�� 0 origname origName
�� 
file�� 0 original  �� 0 repeatalias repeatAlias
�� 
leng�� .0 igorextensionsfldrapp IgorExtensionsFldrApp�� 0 duplicatexops duplicateXOPs
�� .sysobeepnull��� ��� long��o��lv��lvlvE�O���lv��lvlv%E�O*�k+ 
 ��%�%��kva ka ka  OhY hOa �a kva ka ka  Oa r *a a a 0a ,E` W X  *a a l  E` OhO �a !a "a #l $a %%E` &O*a '_ &a (0E` )Oa *E` +OeE` ,O*a 'a !a "a #l $a -%a (0E` .O_ .a 'a //j 0 I*a 1a 'a 2_ .a 3a 4a 5a 6a a 7a  8Oa *a 9,a 'a :/a 6,FO)_ .k+ ;Y hW RX   2_ a #&a <%E` &O*a '_ &a (0E` )Oa =E` +OfE` ,W X  *a a >l  E` )OPOPO_ )a ?,a ?,E` @O*a A_ @a #&a B%/E` CO*a '_ @a #&a D%a (0E` EUO *a A_ )a #&a F%/E` GW !X  a H�a Ikva ka ka  OhO_ +�%a J%�%a K%�%_ a #&%E` +O*_ C_ Gl+ L da  D "*a 1a Ma 2_ Ca N_ Ga  8E` OW !X  a P�a Qkva ka ka  OhUO_ +�%a R%�%a S%E` +OPY _ +�%a T%�%a U%�%a V%E` +OPOa WE` XOa YE` ZO:�[a 1a [l \kh  �a [k/E` ]O�a [l/E` ^O 7_ a #&a _%_ ]%a `%_ ^%E` aOa  *a b_ aa (0E` cUW OX  a d_ ^%a e%�a fkva ka ka  Oa g�%a h%_ a%a i%�a jkva ka ka  OhO*_ E_ cl+ L la  b_ Z�%_ ]%a k%_ ^%E` ZO "*a 1a Ma 2_ Ea N_ ca  8E` lW )X  a m_ ^%a n%�a okva ka ka  OhUY _ X�%_ ^%E` X[OY��O_ Xa p,j _ +�%a q%�%a r%_ X%E` +Y hO_ Za p,j _ +�%a s%�%a t%_ Z%E` +Y hO_ ,e  �a  *a '_ a #&a u%a (0E` vUOa wE` xO ��[a 1a [l \kh  �a [k/E` ]O�a [l/E` ^O_ a #&a y%_ ]%a z%_ ^%E` aOa  *a b_ aa (0E` cUO*_ v_ cl+ L _ x�%_ ^%E` xY h[OY��O_ xa p,j >*j {Oa |_ x%�%�%a }%_ va #&%�%�%a ~%�a kva ka la  Y hY hO_ +�a �kva ka ka  Oh� ����    ����      ����     ! ����    0 3 ����    9 <�  �� �� w��
�� 
sdsk
�� 
cfol �  A p p l i c a t i o n s
�� 
cfol � & I g o r   P r o   6 . 2   F o l d e r� �   � M a c i n t o s h   H D : U s e r s : t i s c h l e r : D o c u m e n t s : W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s : I g o r   S h a r e d   P r o c e d u r e s :� !! "��#" $��%$ &��'& (��)( *��+* ,��-, .��/. w��
�� 
sdsk
�� 
cfol/ �00 
 U s e r s
�� 
cfol- �11  t i s c h l e r
�� 
cfol+ �22  D o c u m e n t s
�� 
cfol) �33  W a v e M e t r i c s
�� 
cfol' �44 * I g o r   P r o   6   U s e r   F i l e s
�� 
cfol% �55  U s e r   P r o c e d u r e s
�� 
cfol# �66 , I g o r   S h a r e d   P r o c e d u r e s� �77 i n s t a l l i n g   i n t o   D o c u m e n t s   F o l d e r  - - - - - - - - - - - - -  u s i n g   I g o r P r o F o l d e r :  M a c i n t o s h   H D : A p p l i c a t i o n s : I g o r   P r o   6 . 2   F o l d e r :  - - - - - - - - - - - - -  ' a l w a y s   a l i a s '   i s   O K ,  n o t h i n g   c h a n g e d .  - - - - - - - - - - - - -  D i d   n o t   c h a n g e   a l i a s e s   t o :  H D F 5   H e l p . i h f  H D F 5 . x o p  H F S A n d P o s i x   H e l p . i h f  H F S A n d P o s i x . x o p
�L boovtrue� 88 9��:9 ;��<; =��>= ?��@? A��BA C��DC w��
�� 
sdsk
�� 
cfolD �EE 
 U s e r s
�� 
cfolB �FF  t i s c h l e r
�� 
cfol@ �GG  D o c u m e n t s
�� 
cfol> �HH  W a v e M e t r i c s
�� 
cfol< �II * I g o r   P r o   6   U s e r   F i l e s
�� 
cfol: �JJ  U s e r   P r o c e d u r e s� KK L��ML N��ON P��QP R��SR T��UT w��
�� 
sdsk
�� 
cfolU �VV 
 U s e r s
�� 
cfolS �WW  t i s c h l e r
�� 
cfolQ �XX  D o c u m e n t s
�� 
cfolO �YY  W a v e M e t r i c s
�� 
cfolM �ZZ * I g o r   P r o   6   U s e r   F i l e s� [[ \��]\ ^��_^ `��a` b��cb d��ed f��gf w��
�� 
sdsk
�� 
cfolg �hh 
 U s e r s
�� 
cfole �ii  t i s c h l e r
�� 
cfolc �jj  D o c u m e n t s
�� 
cfola �kk  W a v e M e t r i c s
�� 
cfol_ �ll * I g o r   P r o   6   U s e r   F i l e s
�� 
cfol] �mm  I g o r   P r o c e d u r e s� nn o��po q��rq s��ts u��vu w��xw y��zy w��
�� 
sdsk
�� 
cfolz �{{ 
 U s e r s
�� 
cfolx �||  t i s c h l e r
�� 
cfolv �}}  D o c u m e n t s
�� 
cfolt �~~  W a v e M e t r i c s
�� 
cfolr � * I g o r   P r o   6   U s e r   F i l e s
�� 
cfolp ���  I g o r   E x t e n s i o n s�alis      Macintosh HD               Ũ�H+   ��always                                                          ��ť��        ����  	                Igor Shared Procedures    ũz      ŦJ     �� B� B� B� 1� 	I�  {}  uMacintosh HD:Users:tischler:Documents:WaveMetrics:Igor Pro 6 User Files:User Procedures:Igor Shared Procedures:always     a l w a y s    M a c i n t o s h   H D  hUsers/tischler/Documents/WaveMetrics/Igor Pro 6 User Files/User Procedures/Igor Shared Procedures/always  /    ��  � ��� x  H D F 5   H e l p . i h f  H D F 5 . x o p  H F S A n d P o s i x   H e l p . i h f  H F S A n d P o s i x . x o p� ��� � M a c i n t o s h   H D : A p p l i c a t i o n s : I g o r   P r o   6 . 2   F o l d e r : M o r e   E x t e n s i o n s : U t i l i t i e s : H F S A n d P o s i x . x o p� �� ����� ����� ����� ����� ����� w��
�� 
sdsk
�� 
cfol� ���  A p p l i c a t i o n s
�� 
cfol� ��� & I g o r   P r o   6 . 2   F o l d e r
�� 
cfol� ���  M o r e   E x t e n s i o n s
�� 
cfol� ���  U t i l i t i e s
�� 
docf� ���  H F S A n d P o s i x . x o p� �� ���� ��~�� ��}�� w�|
�| 
sdsk
�} 
cfol� ���  A p p l i c a t i o n s
�~ 
cfol� ��� & I g o r   P r o   6 . 2   F o l d e r
� 
cfol� ���  I g o r   E x t e n s i o n s�K  �J  �I  �H  �G  �F  �E  �D  �C  ascr  ��ޭ