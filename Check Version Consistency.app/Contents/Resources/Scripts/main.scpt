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
currentDir��    m      " "�                                                                                  MACS  alis    n  
JZTdesktop                 ΍��H+     {
Finder.app                                                      %��_��        ����  	                CoreServices    ΍�      �`D       {   x   w  4JZTdesktop:System: Library: CoreServices: Finder.app   
 F i n d e r . a p p   
 J Z T d e s k t o p  &System/Library/CoreServices/Finder.app  / ��  ��  ��     # $ # l    %���� % r     & ' & l    (���� ( c     ) * ) l    +���� + I   �� ,��
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
currentDir��  ��   V m   G J k k � l l  " T o      ���� 0 posixcmd PosixCmd��  ��   Q  m n m l  P _ o���� o I  P _�� p��
�� .sysodlogaskr        TEXT p b   P [ q r q b   P W s t s m   P S u u � v v  P O S I X   p a t h : t o   S V��
�� 
ret  r l  W Z w���� w n   W Z x y x 1   X Z��
�� 
psxp y o   W X���� 0 
currentdir 
currentDir��  ��  ��  ��  ��   n  z { z l     �� | }��   |  log PosixCmd    } � ~ ~  l o g   P o s i x C m d {   �  l     ��������  ��  ��   �  � � � l  ` k ����� � r   ` k � � � I  ` g�� ���
�� .sysoexecTEXT���     TEXT � o   ` c���� 0 posixcmd PosixCmd��   � o      ���� 0 vstatus vStatus��  ��   �  � � � l     �� � ���   �  log vStatus    � � � �  l o g   v S t a t u s �  � � � l     ��������  ��  ��   �  � � � l     ��������  ��  ��   �  � � � l  l � ����� � Z   l � � ��� � � C   l s � � � o   l o���� 0 vstatus vStatus � m   o r � � � � �  P e r f e c t   M a t c h � k   v � � �  � � � l  v v�� � ���   � N Hdisplay dialog vStatus buttons {"Continue"} default button 1 with icon 1    � � � � � d i s p l a y   d i a l o g   v S t a t u s   b u t t o n s   { " C o n t i n u e " }   d e f a u l t   b u t t o n   1   w i t h   i c o n   1 �  ��� � I  v ��� � �
�� .sysodlogaskr        TEXT � o   v y���� 0 vstatus vStatus � �� � �
�� 
btns � J   | � � �  ��� � m   |  � � � � �  C o n t i n u e��   � �� ���
�� 
dflt � m   � ����� ��  ��  ��   � O   � � � � � k   � � � �  � � � I  � �������
�� .miscactvnull��� ��� null��  ��   �  � � � Z   � � � ��� � � =   � � � � � l  � � ����� � I  � ��� ���
�� .coredoexbool       obj  � 4  � ��� �
�� 
docu � m   � ����� ��  ��  ��   � m   � ���
�� boovtrue � Z   � � � ����� � >   � � � � � l  � � ����� � n   � � � � � m   � ���
�� 
ctxt � l  � � ����� � 4  � ��� �
�� 
docu � m   � ����� ��  ��  ��  ��   � m   � � � � � � �   � I  � ����� �
�� .corecrel****      � null��   � �� ��
�� 
kocl � m   � ��~
�~ 
docu�  ��  ��  ��   � I  � ��}�| �
�} .corecrel****      � null�|   � �{ ��z
�{ 
kocl � m   � ��y
�y 
docu�z   �  ��x � r   � � � � � o   � ��w�w 0 vstatus vStatus � l      ��v�u � n       � � � m   � ��t
�t 
ctxt � l  � � ��s�r � 4  � ��q �
�q 
docu � m   � ��p�p �s  �r  �v  �u  �x   � 4   � ��o �
�o 
capp � I   � ��n�m�l�n 00 displayapplicationname displayApplicationName�m  �l  ��  ��   �  � � � l     �k�j�i�k  �j  �i   �  � � � l     �h�g�f�h  �g  �f   �  � � � l     �e � ��e   � C = return name of application to use for displaying mis-matches    � � � � z   r e t u r n   n a m e   o f   a p p l i c a t i o n   t o   u s e   f o r   d i s p l a y i n g   m i s - m a t c h e s �  � � � i      � � � I      �d�c�b�d 00 displayapplicationname displayApplicationName�c  �b   � k     6 � �  � � � l     �a � ��a   � 8 2 prefered app is BBEdit, if it is on this computer    � � � � d   p r e f e r e d   a p p   i s   B B E d i t ,   i f   i t   i s   o n   t h i s   c o m p u t e r �  � � � r      � � � m      � � � � �  B B E d i t � o      �`�` 0 
displayapp 
displayApp �  � � � Q     � � � � I   �_ ��^
�_ .sysoexecTEXT���     TEXT � b     � � � b    
 � � � m     � � � � � D o s a s c r i p t   - e   ' e x i s t s   a p p l i c a t i o n   " � o    	�]�] 0 
displayapp 
displayApp � m   
  � � � � �  " '�^   � R      �\�[ �
�\ .ascrerr ****      � ****�[   � �Z ��Y
�Z 
errn � o      �X�X 0 	errnumber 	errNumber�Y   � k     � �  � � � l   �W � ��W   � ) # BBEdit not found, try TextWrangler    � � � � F   B B E d i t   n o t   f o u n d ,   t r y   T e x t W r a n g l e r �  �V  r     m     �  T e x t W r a n g l e r o      �U�U 0 
displayapp 
displayApp�V   �  l   �T�S�R�T  �S  �R    Q    3	
	 I   (�Q�P
�Q .sysoexecTEXT���     TEXT b    $ b    " m      � D o s a s c r i p t   - e   ' e x i s t s   a p p l i c a t i o n   " o     !�O�O 0 
displayapp 
displayApp m   " # �  " '�P  
 R      �N�M
�N .ascrerr ****      � ****�M   �L�K
�L 
errn o      �J�J 0 	errnumber 	errNumber�K   k   0 3  l  0 0�I�I   P J TextWrangler failed too, just use TextEdit, it should always be available    � �   T e x t W r a n g l e r   f a i l e d   t o o ,   j u s t   u s e   T e x t E d i t ,   i t   s h o u l d   a l w a y s   b e   a v a i l a b l e �H r   0 3 m   0 1   �!!  T e x t E d i t o      �G�G 0 
displayapp 
displayApp�H   "�F" L   4 6## o   4 5�E�E 0 
displayapp 
displayApp�F   � $�D$ l     �C�B�A�C  �B  �A  �D       
�@%&'()*+,�?�@  % �>�=�<�;�:�9�8�7�> 00 displayapplicationname displayApplicationName
�= .aevtoappnull  �   � ****�< 0 
currentdir 
currentDir�;  0 applescriptapp AppleScriptApp�: 0 	pythonmac 	pythonMac�9 0 posixcmd PosixCmd�8 0 vstatus vStatus�7  & �6 ��5�4-.�3�6 00 displayapplicationname displayApplicationName�5  �4  - �2�1�2 0 
displayapp 
displayApp�1 0 	errnumber 	errNumber. 
 � � ��0�// 
�0 .sysoexecTEXT���     TEXT�/  / �.�-�,
�. 
errn�- 0 	errnumber 	errNumber�,  �3 7�E�O �%�%j W 
X  �E�O �%�%j W 
X  �E�O�' �+0�*�)12�(
�+ .aevtoappnull  �   � ****0 k     �33  44  #55  -66  ?77  P88  m99  �::  ��'�'  �*  �)  1  2 $ "�&�%�$�#�" 4�!� �� K� _� d f k� u��� �� �������� ���
�& .earsffdralis        afdr
�% 
ctnr
�$ 
ctxt�# 0 
currentdir 
currentDir�"  0 applescriptapp AppleScriptApp�! 0 	pythonmac 	pythonMac
�  
alis�  �  
� .sysodlogaskr        TEXT
� 
psxp� 0 posixcmd PosixCmd
� 
ret 
� .sysoexecTEXT���     TEXT� 0 vstatus vStatus
� 
btns
� 
dflt� 
� 
capp� 00 displayapplicationname displayApplicationName
� .miscactvnull��� ��� null
� 
docu
� .coredoexbool       obj 
� 
kocl
� .corecrel****      � null�( �� )j �,�&E�UO)j �&E�O��%E�O ��&W X 	 
�j OhO���,%�%a %��,%a %E` Oa _ %��,%j O_ j E` O_ a  _ a a kva ka  Y ^*a *j+ / P*j O*a k/j  e  #*a k/�-a ! *a "a l #Y hY *a "a l #O_ *a k/�-FU( �;; � t i s c h l e r : D o c u m e n t s : W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s : L o c a l P a c k a g e s :) �<< � t i s c h l e r : D o c u m e n t s : W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s : L o c a l P a c k a g e s : C h e c k   V e r s i o n   C o n s i s t e n c y . a p p :* �==2 t i s c h l e r : D o c u m e n t s : W a v e M e t r i c s : I g o r   P r o   6   U s e r   F i l e s : U s e r   P r o c e d u r e s : L o c a l P a c k a g e s : C h e c k   V e r s i o n   C o n s i s t e n c y . a p p : C o n t e n t s : R e s o u r c e s : c h e c k V e r s i o n S t a t u s . p y+ �>> " / V o l u m e s / t i s c h l e r / D o c u m e n t s / W a v e M e t r i c s / I g o r   P r o   6   U s e r   F i l e s / U s e r   P r o c e d u r e s / L o c a l P a c k a g e s / C h e c k   V e r s i o n   C o n s i s t e n c y . a p p / C o n t e n t s / R e s o u r c e s / c h e c k V e r s i o n S t a t u s . p y "   " / V o l u m e s / t i s c h l e r / D o c u m e n t s / W a v e M e t r i c s / I g o r   P r o   6   U s e r   F i l e s / U s e r   P r o c e d u r e s / L o c a l P a c k a g e s / ", �?? . P e r f e c t   M a t c h ,  6 5   f i l e s�?   ascr  ��ޭ