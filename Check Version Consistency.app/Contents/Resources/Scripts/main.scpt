FasdUAS 1.101.10   ��   ��    k             l      ��  ��    �  Check the consitency of VersionStatus.xml with the ipf files acutally present
This just is an easy way to run the python file
     � 	 	 �   C h e c k   t h e   c o n s i t e n c y   o f   V e r s i o n S t a t u s . x m l   w i t h   t h e   i p f   f i l e s   a c u t a l l y   p r e s e n t 
 T h i s   j u s t   i s   a n   e a s y   w a y   t o   r u n   t h e   p y t h o n   f i l e 
   
  
 l     ��������  ��  ��        l     ����  O         k           l   ��  ��    N Hset currentDir to the current working directory, where this is run from.     �   � s e t   c u r r e n t D i r   t o   t h e   c u r r e n t   w o r k i n g   d i r e c t o r y ,   w h e r e   t h i s   i s   r u n   f r o m .   ��  r        l    ����  c        l    ����  n        m   	 ��
�� 
ctnr  l   	  ����   I   	�� !��
�� .earsffdralis        afdr !  f    ��  ��  ��  ��  ��    m    ��
�� 
ctxt��  ��    o      ���� 0 
currentdir 
currentDir��    m      " "�                                                                                  MACS  alis    8  Spurr HD                       BD ����
Finder.app                                                     ����            ����  
 cu             CoreServices  )/:System:Library:CoreServices:Finder.app/    
 F i n d e r . a p p    S p u r r   H D  &System/Library/CoreServices/Finder.app  / ��  ��  ��     # $ # l    %���� % r     & ' & l    (���� ( c     ) * ) l    +���� + I   �� ,��
�� .earsffdralis        afdr ,  f    ��  ��  ��   * m    ��
�� 
ctxt��  ��   ' o      ����  0 applescriptapp AppleScriptApp��  ��   $  - . - l     /���� / r      0 1 0 b     2 3 2 o    ����  0 applescriptapp AppleScriptApp 3 m     4 4 � 5 5 P C o n t e n t s : R e s o u r c e s : c h e c k V e r s i o n S t a t u s . p y 1 o      ���� 0 	pythonmac 	pythonMac��  ��   .  6 7 6 l     ��������  ��  ��   7  8 9 8 l     ��������  ��  ��   9  : ; : l     �� < =��   < 2 , check on existance of checkVersionStatus.py    = � > > X   c h e c k   o n   e x i s t a n c e   o f   c h e c k V e r s i o n S t a t u s . p y ;  ? @ ? l  ! 7 A���� A Q   ! 7 B C D B c   $ ' E F E o   $ %���� 0 	pythonmac 	pythonMac F m   % &��
�� 
alis C R      ������
�� .ascrerr ****      � ****��  ��   D k   / 7 G G  H I H I  / 4�� J��
�� .sysodlogaskr        TEXT J m   / 0 K K � L L B C a n n o t   f i n d   c h e c k V e r s i o n S t a t u s . p y��   I  M�� M L   5 7����  ��  ��  ��   @  N O N l     ��������  ��  ��   O  P Q P l  8 O R���� R r   8 O S T S b   8 K U V U b   8 G W X W b   8 C Y Z Y b   8 ? [ \ [ b   8 = ] ^ ] m   8 9 _ _ � ` `  " ^ l  9 < a���� a n   9 < b c b 1   : <��
�� 
psxp c o   9 :���� 0 	pythonmac 	pythonMac��  ��   \ m   = > d d � e e  " Z m   ? B f f � g g    " X l  C F h���� h n   C F i j i 1   D F��
�� 
psxp j o   C D���� 0 
currentdir 
currentDir��  ��   V m   G J k k � l l  " T o      ���� 0 posixcmd PosixCmd��  ��   Q  m n m l     �� o p��   o  log PosixCmd    p � q q  l o g   P o s i x C m d n  r s r l     ��������  ��  ��   s  t u t l  P [ v���� v r   P [ w x w I  P W�� y��
�� .sysoexecTEXT���     TEXT y o   P S���� 0 posixcmd PosixCmd��   x o      ���� 0 vstatus vStatus��  ��   u  z { z l     �� | }��   |  log vStatus    } � ~ ~  l o g   v S t a t u s {   �  l     ��������  ��  ��   �  � � � l     ��������  ��  ��   �  � � � l  \ � ����� � Z   \ � � ��� � � C   \ c � � � o   \ _���� 0 vstatus vStatus � m   _ b � � � � �  P e r f e c t   M a t c h � k   f { � �  � � � l  f f�� � ���   � N Hdisplay dialog vStatus buttons {"Continue"} default button 1 with icon 1    � � � � � d i s p l a y   d i a l o g   v S t a t u s   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1 �  ��� � I  f {�� � �
�� .sysodlogaskr        TEXT � o   f i���� 0 vstatus vStatus � �� � �
�� 
btns � J   l q � �  ��� � m   l o � � � � �  C o n t i n u e��   � �� ���
�� 
dflt � m   t u���� ��  ��  ��   � O   ~ � � � � k   � � � �  � � � I  � �������
�� .miscactvnull��� ��� null��  ��   �  � � � Z   � � � ��� � � =   � � � � � l  � � ����� � I  � ��� ���
�� .coredoexbool       obj  � 4  � ��� �
�� 
docu � m   � ����� ��  ��  ��   � m   � ���
�� boovtrue � Z   � � � ����� � >   � � � � � l  � � ����� � n   � � � � � m   � ���
�� 
ctxt � l  � � ����� � 4  � ��� �
�� 
docu � m   � ����� ��  ��  ��  ��   � m   � � � � � � �   � I  � ����� �
�� .corecrel****      � null��   � �� ���
�� 
kocl � m   � ���
�� 
docu��  ��  ��  ��   � I  � ����� �
�� .corecrel****      � null��   � �� ���
�� 
kocl � m   � ���
�� 
docu��   �  ��� � r   � � � � � o   � ����� 0 vstatus vStatus � l      ���~ � n       � � � m   � ��}
�} 
ctxt � l  � � ��|�{ � 4  � ��z �
�z 
docu � m   � ��y�y �|  �{  �  �~  ��   � 4   ~ ��x �
�x 
capp � I   � ��w�v�u�w 00 displayapplicationname displayApplicationName�v  �u  ��  ��   �  � � � l     �t�s�r�t  �s  �r   �  � � � l     �q�p�o�q  �p  �o   �  � � � l     �n � ��n   � C = return name of application to use for displaying mis-matches    � � � � z   r e t u r n   n a m e   o f   a p p l i c a t i o n   t o   u s e   f o r   d i s p l a y i n g   m i s - m a t c h e s �  � � � i      � � � I      �m�l�k�m 00 displayapplicationname displayApplicationName�l  �k   � k     6 � �  � � � l     �j � ��j   � 8 2 prefered app is BBEdit, if it is on this computer    � � � � d   p r e f e r e d   a p p   i s   B B E d i t ,   i f   i t   i s   o n   t h i s   c o m p u t e r �  � � � r      � � � m      � � � � �  B B E d i t � o      �i�i 0 
displayapp 
displayApp �  � � � Q     � � � � I   �h ��g
�h .sysoexecTEXT���     TEXT � b     � � � b    
 � � � m     � � � � � D o s a s c r i p t   - e   ' e x i s t s   a p p l i c a t i o n   " � o    	�f�f 0 
displayapp 
displayApp � m   
  � � � � �  " '�g   � R      �e�d �
�e .ascrerr ****      � ****�d   � �c ��b
�c 
errn � o      �a�a 0 	errnumber 	errNumber�b   � k     � �  � � � l   �` � ��`   � ) # BBEdit not found, try TextWrangler    � � � � F   B B E d i t   n o t   f o u n d ,   t r y   T e x t W r a n g l e r �  ��_ � r     � � � m     � � � � �  T e x t W r a n g l e r � o      �^�^ 0 
displayapp 
displayApp�_   �  � � � l   �]�\�[�]  �\  �[   �  � � � Q    3 � � � � I   (�Z ��Y
�Z .sysoexecTEXT���     TEXT � b    $   b    " m      � D o s a s c r i p t   - e   ' e x i s t s   a p p l i c a t i o n   " o     !�X�X 0 
displayapp 
displayApp m   " # �  " '�Y   � R      �W�V
�W .ascrerr ****      � ****�V   �U	�T
�U 
errn	 o      �S�S 0 	errnumber 	errNumber�T   � k   0 3

  l  0 0�R�R   P J TextWrangler failed too, just use TextEdit, it should always be available    � �   T e x t W r a n g l e r   f a i l e d   t o o ,   j u s t   u s e   T e x t E d i t ,   i t   s h o u l d   a l w a y s   b e   a v a i l a b l e �Q r   0 3 m   0 1 �  T e x t E d i t o      �P�P 0 
displayapp 
displayApp�Q   � �O L   4 6 o   4 5�N�N 0 
displayapp 
displayApp�O   � �M l     �L�K�J�L  �K  �J  �M       �I�I   �H�G�H 00 displayapplicationname displayApplicationName
�G .aevtoappnull  �   � **** �F ��E�D�C�F 00 displayapplicationname displayApplicationName�E  �D   �B�A�B 0 
displayapp 
displayApp�A 0 	errnumber 	errNumber 
 � � ��@�? �
�@ .sysoexecTEXT���     TEXT�?   �>�=�<
�> 
errn�= 0 	errnumber 	errNumber�<  �C 7�E�O �%�%j W 
X  �E�O �%�%j W 
X  �E�O� �;�:�9 �8
�; .aevtoappnull  �   � **** k     �!!  ""  ###  -$$  ?%%  P&&  t''  ��7�7  �:  �9      " "�6�5�4�3�2 4�1�0�/�. K�- _�, d f k�+�*�) ��( ��'�&�%�$�#�"�! �� �
�6 .earsffdralis        afdr
�5 
ctnr
�4 
ctxt�3 0 
currentdir 
currentDir�2  0 applescriptapp AppleScriptApp�1 0 	pythonmac 	pythonMac
�0 
alis�/  �.  
�- .sysodlogaskr        TEXT
�, 
psxp�+ 0 posixcmd PosixCmd
�* .sysoexecTEXT���     TEXT�) 0 vstatus vStatus
�( 
btns
�' 
dflt�& 
�% 
capp�$ 00 displayapplicationname displayApplicationName
�# .miscactvnull��� ��� null
�" 
docu
�! .coredoexbool       obj 
�  
kocl
� .corecrel****      � null�8 �� )j �,�&E�UO)j �&E�O��%E�O ��&W X 	 
�j OhO���,%�%a %��,%a %E` O_ j E` O_ a  _ a a kva ka  Y ^*a *j+ / P*j O*a k/j e  #*a k/�-a  *a  a l !Y hY *a  a l !O_ *a k/�-FUascr  ��ޭ