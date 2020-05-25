!**************************************************************
!* AceGen    6.702 Windows (4 May 16)                         *
!*           Co. J. Korelc  2013           10 Aug 19 21:43:00 *
!**************************************************************
! User     : Full professional version
! Notebook : MainFile
! Evaluation time                 : 23 s    Mode  : Optimal
! Number of formulae              : 620     Method: Automatic
! Subroutine                      : elastic size: 13814
! Total size of Mathematica  code : 13814 subexpressions
! Total size of Fortran code      : 31084 bytes

!******************* S U B R O U T I N E **********************
SUBROUTINE elastic(v,mpar,statev,Fnew,sigma,ddsdde,yielding,xguess)
USE SMSUtility
IMPLICIT NONE
INTEGER i01,i02
DOUBLE PRECISION v(871),mpar(21),statev(36),Fnew(9),sigma(6),ddsdde(6,6),yielding,xguess(16)
v(858)=mpar(3)+statev(35)+mpar(6)*(1d0-dexp(-(statev(10)/mpar(4))))+mpar(7)*(1d0-dexp(-(statev(10)/mpar(5))))
v(857)=2d0*statev(33)**2
v(856)=2d0*statev(32)**2
v(855)=2d0*statev(31)**2
v(854)=statev(30)**2
v(853)=statev(29)**2
v(825)=-statev(29)-statev(30)
v(824)=statev(7)*statev(9)
v(823)=-(statev(6)*statev(7))
v(822)=statev(5)*statev(6)
v(821)=-(statev(5)*statev(9))
v(820)=1d0+statev(3)
v(819)=statev(4)*statev(5)
v(818)=-(statev(4)*statev(8))
v(817)=statev(8)*statev(9)
v(816)=1d0+statev(2)
v(815)=statev(4)*statev(6)
v(814)=statev(7)*statev(8)
v(813)=1d0+statev(1)
v(812)=statev(24)*statev(25)
v(811)=-(statev(25)*statev(26))
v(810)=statev(26)*statev(28)
v(809)=-(statev(24)*statev(28))
v(808)=1d0+statev(22)
v(807)=statev(27)*statev(28)
v(806)=-(statev(23)*statev(27))
v(805)=statev(23)*statev(24)
v(804)=1d0+statev(21)
v(803)=statev(23)*statev(25)
v(802)=statev(26)*statev(27)
v(801)=1d0+statev(20)
v(800)=statev(15)*statev(16)
v(799)=-(statev(16)*statev(17))
v(798)=statev(17)*statev(19)
v(797)=-(statev(15)*statev(19))
v(796)=1d0+statev(13)
v(795)=statev(18)*statev(19)
v(794)=-(statev(14)*statev(18))
v(793)=statev(14)*statev(15)
v(792)=1d0+statev(12)
v(791)=statev(14)*statev(16)
v(790)=statev(17)*statev(18)
v(789)=1d0+statev(11)
v(788)=1d0/(Fnew(6)*(Fnew(4)*Fnew(5)-Fnew(2)*Fnew(7))+Fnew(3)*(Fnew(1)*Fnew(2)-Fnew(4)*Fnew(8))+(-(Fnew(1)*Fnew(5))&
&+Fnew(7)*Fnew(8))*Fnew(9))
v(863)=2d0*v(788)
v(180)=-(statev(15)*v(789))+v(790)
v(171)=-(statev(19)*v(789))+v(791)
v(179)=-(statev(17)*v(792))+v(793)
v(178)=v(789)*v(792)+v(794)
v(172)=-(statev(16)*v(792))+v(795)
v(176)=v(792)*v(796)+v(797)
v(175)=-(statev(14)*v(796))+v(798)
v(174)=v(789)*v(796)+v(799)
v(173)=-(statev(18)*v(796))+v(800)
v(203)=-(statev(24)*v(801))+v(802)
v(194)=-(statev(28)*v(801))+v(803)
v(202)=-(statev(26)*v(804))+v(805)
v(201)=v(801)*v(804)+v(806)
v(195)=-(statev(25)*v(804))+v(807)
v(199)=v(804)*v(808)+v(809)
v(198)=-(statev(23)*v(808))+v(810)
v(197)=v(801)*v(808)+v(811)
v(196)=-(statev(27)*v(808))+v(812)
v(85)=-(statev(5)*v(813))+v(814)
v(81)=-(statev(9)*v(813))+v(815)
v(89)=-(statev(6)*v(816))+v(817)
v(87)=v(813)*v(816)+v(818)
v(86)=-(statev(7)*v(816))+v(819)
v(91)=v(816)*v(820)+v(821)
v(90)=-(statev(8)*v(820))+v(822)
v(83)=v(813)*v(820)+v(823)
v(82)=-(statev(4)*v(820))+v(824)
v(77)=1d0/(statev(9)*v(85)+statev(6)*v(86)+v(820)*v(87))
v(240)=v(77)*v(85)
v(239)=v(77)*v(86)
v(238)=v(77)*v(87)
v(237)=v(77)*v(82)
v(236)=v(77)*v(81)
v(235)=v(77)*v(83)
v(234)=v(77)*v(89)
v(233)=v(77)*v(90)
v(232)=v(77)*v(91)
v(78)=v(77)*(Fnew(7)*v(89)+Fnew(4)*v(90)+Fnew(1)*v(91))
v(859)=(v(78)*v(78))
v(826)=2d0*v(78)
v(247)=v(234)*v(826)
v(256)=-v(247)/3d0
v(244)=v(233)*v(826)
v(253)=-v(244)/3d0
v(241)=v(232)*v(826)
v(250)=-v(241)/3d0
v(79)=v(77)*(Fnew(5)*v(81)+Fnew(8)*v(82)+Fnew(2)*v(83))
v(827)=2d0*v(79)
v(266)=v(237)*v(827)
v(275)=-v(266)/3d0
v(263)=v(236)*v(827)
v(272)=-v(263)/3d0
v(260)=v(235)*v(827)
v(269)=-v(260)/3d0
v(80)=v(77)*(Fnew(9)*v(85)+Fnew(6)*v(86)+Fnew(3)*v(87))
v(828)=2d0*v(80)
v(285)=v(240)*v(828)
v(294)=-v(285)/3d0
v(282)=v(239)*v(828)
v(291)=-v(282)/3d0
v(279)=v(238)*v(828)
v(288)=-v(279)/3d0
v(84)=v(77)*(Fnew(7)*v(81)+Fnew(1)*v(82)+Fnew(4)*v(83))
v(860)=(v(84)*v(84))
v(829)=2d0*v(84)
v(301)=v(236)*v(78)+v(234)*v(84)
v(298)=v(235)*v(78)+v(233)*v(84)
v(295)=v(237)*v(78)+v(232)*v(84)
v(265)=v(236)*v(829)
v(274)=-v(265)/3d0
v(262)=v(235)*v(829)
v(271)=-v(262)/3d0
v(259)=v(237)*v(829)
v(268)=-v(259)/3d0
v(88)=v(77)*(Fnew(2)*v(85)+Fnew(8)*v(86)+Fnew(5)*v(87))
v(830)=2d0*v(88)
v(329)=v(239)*v(79)+v(237)*v(88)
v(326)=v(238)*v(79)+v(236)*v(88)
v(323)=v(240)*v(79)+v(235)*v(88)
v(284)=v(239)*v(830)
v(293)=-v(284)/3d0
v(281)=v(238)*v(830)
v(290)=-v(281)/3d0
v(278)=v(240)*v(830)
v(287)=-v(278)/3d0
v(92)=v(77)*(Fnew(3)*v(89)+Fnew(9)*v(90)+Fnew(6)*v(91))
v(831)=2d0*v(92)
v(357)=v(233)*v(80)+v(240)*v(92)
v(354)=v(232)*v(80)+v(239)*v(92)
v(351)=v(234)*v(80)+v(238)*v(92)
v(249)=v(233)*v(831)
v(258)=-v(249)/3d0
v(246)=v(232)*v(831)
v(255)=-v(246)/3d0
v(243)=v(234)*v(831)
v(252)=-v(243)/3d0
v(93)=v(77)*(Fnew(4)*v(85)+Fnew(1)*v(86)+Fnew(7)*v(87))
v(862)=v(829)*v(93)
v(861)=(v(93)*v(93))
v(832)=2d0*v(93)
v(355)=v(238)*v(78)+v(234)*v(93)
v(352)=v(240)*v(78)+v(233)*v(93)
v(349)=v(239)*v(78)+v(232)*v(93)
v(328)=v(238)*v(84)+v(236)*v(93)
v(325)=v(240)*v(84)+v(235)*v(93)
v(322)=v(239)*v(84)+v(237)*v(93)
v(283)=v(238)*v(832)
v(292)=-v(283)/3d0
v(280)=v(240)*v(832)
v(289)=-v(280)/3d0
v(277)=v(239)*v(832)
v(286)=-v(277)/3d0
v(94)=v(77)*(Fnew(5)*v(89)+Fnew(2)*v(90)+Fnew(8)*v(91))
v(833)=2d0*v(94)
v(356)=v(232)*v(88)+v(239)*v(94)
v(353)=v(234)*v(88)+v(238)*v(94)
v(350)=v(233)*v(88)+v(240)*v(94)
v(302)=v(232)*v(79)+v(237)*v(94)
v(299)=v(234)*v(79)+v(236)*v(94)
v(296)=v(233)*v(79)+v(235)*v(94)
v(248)=v(232)*v(833)
v(257)=-v(248)/3d0
v(245)=v(234)*v(833)
v(254)=-v(245)/3d0
v(242)=v(233)*v(833)
v(251)=-v(242)/3d0
v(95)=v(77)*(Fnew(3)*v(81)+Fnew(6)*v(82)+Fnew(9)*v(83))
v(834)=2d0*v(95)
v(330)=v(235)*v(80)+v(240)*v(95)
v(327)=v(237)*v(80)+v(239)*v(95)
v(324)=v(236)*v(80)+v(238)*v(95)
v(303)=v(235)*v(92)+v(233)*v(95)
v(300)=v(237)*v(92)+v(232)*v(95)
v(297)=v(236)*v(92)+v(234)*v(95)
v(267)=v(235)*v(834)
v(276)=-v(267)/3d0
v(264)=v(237)*v(834)
v(273)=-v(264)/3d0
v(261)=v(236)*v(834)
v(270)=-v(261)/3d0
v(96)=v(859)+(v(92)*v(92))+(v(94)*v(94))
v(116)=-v(96)/3d0
v(97)=(v(79)*v(79))+v(860)+(v(95)*v(95))
v(117)=-v(97)/3d0
v(98)=(v(80)*v(80))+v(861)+(v(88)*v(88))
v(446)=v(116)+v(117)+(2d0/3d0)*v(98)
v(109)=-v(98)/3d0
v(456)=v(109)+v(116)+(2d0/3d0)*v(97)
v(436)=v(109)+v(117)+(2d0/3d0)*v(96)
v(99)=v(78)*v(84)+v(79)*v(94)+v(92)*v(95)
v(835)=2d0*v(99)
v(312)=v(303)*v(835)
v(321)=-v(312)+v(267)*v(96)+v(249)*v(97)
v(311)=v(302)*v(835)
v(320)=-v(311)+v(266)*v(96)+v(248)*v(97)
v(310)=v(301)*v(835)
v(319)=-v(310)+v(265)*v(96)+v(247)*v(97)
v(309)=v(300)*v(835)
v(318)=-v(309)+v(264)*v(96)+v(246)*v(97)
v(308)=v(299)*v(835)
v(317)=-v(308)+v(263)*v(96)+v(245)*v(97)
v(307)=v(298)*v(835)
v(316)=-v(307)+v(262)*v(96)+v(244)*v(97)
v(306)=v(297)*v(835)
v(315)=-v(306)+v(261)*v(96)+v(243)*v(97)
v(305)=v(296)*v(835)
v(314)=-v(305)+v(260)*v(96)+v(242)*v(97)
v(304)=v(295)*v(835)
v(313)=-v(304)+v(259)*v(96)+v(241)*v(97)
v(115)=(v(99)*v(99))
v(131)=-v(115)+v(96)*v(97)
v(100)=v(79)*v(88)+v(84)*v(93)+v(80)*v(95)
v(836)=2d0*v(100)
v(406)=v(100)*v(835)
v(339)=v(330)*v(836)
v(348)=-v(339)+v(285)*v(97)+v(267)*v(98)
v(338)=v(329)*v(836)
v(347)=-v(338)+v(284)*v(97)+v(266)*v(98)
v(337)=v(328)*v(836)
v(346)=-v(337)+v(283)*v(97)+v(265)*v(98)
v(336)=v(327)*v(836)
v(345)=-v(336)+v(282)*v(97)+v(264)*v(98)
v(335)=v(326)*v(836)
v(344)=-v(335)+v(281)*v(97)+v(263)*v(98)
v(334)=v(325)*v(836)
v(343)=-v(334)+v(280)*v(97)+v(262)*v(98)
v(333)=v(324)*v(836)
v(342)=-v(333)+v(279)*v(97)+v(261)*v(98)
v(332)=v(323)*v(836)
v(341)=-v(332)+v(278)*v(97)+v(260)*v(98)
v(331)=v(322)*v(836)
v(340)=-v(331)+v(277)*v(97)+v(259)*v(98)
v(103)=(v(100)*v(100))
v(120)=-v(103)+v(97)*v(98)
v(101)=v(80)*v(92)+v(78)*v(93)+v(88)*v(94)
v(837)=2d0*v(101)
v(405)=v(101)*v(835)
v(403)=v(101)*v(836)
v(393)=v(357)*v(837)
v(402)=-v(393)+v(285)*v(96)+v(249)*v(98)
v(392)=v(356)*v(837)
v(401)=-v(392)+v(284)*v(96)+v(248)*v(98)
v(391)=v(355)*v(837)
v(400)=-v(391)+v(283)*v(96)+v(247)*v(98)
v(390)=v(354)*v(837)
v(399)=-v(390)+v(282)*v(96)+v(246)*v(98)
v(389)=v(353)*v(837)
v(398)=-v(389)+v(281)*v(96)+v(245)*v(98)
v(388)=v(352)*v(837)
v(397)=-v(388)+v(280)*v(96)+v(244)*v(98)
v(387)=v(351)*v(837)
v(396)=-v(387)+v(279)*v(96)+v(243)*v(98)
v(386)=v(350)*v(837)
v(395)=-v(386)+v(278)*v(96)+v(242)*v(98)
v(385)=v(349)*v(837)
v(394)=-v(385)+v(277)*v(96)+v(241)*v(98)
v(384)=v(101)*v(330)+v(100)*v(357)-v(303)*v(98)-v(285)*v(99)
v(383)=v(101)*v(329)+v(100)*v(356)-v(302)*v(98)-v(284)*v(99)
v(382)=v(101)*v(328)+v(100)*v(355)-v(301)*v(98)-v(283)*v(99)
v(381)=v(101)*v(327)+v(100)*v(354)-v(300)*v(98)-v(282)*v(99)
v(380)=v(101)*v(326)+v(100)*v(353)-v(299)*v(98)-v(281)*v(99)
v(379)=v(101)*v(325)+v(100)*v(352)-v(298)*v(98)-v(280)*v(99)
v(378)=v(101)*v(324)+v(100)*v(351)-v(297)*v(98)-v(279)*v(99)
v(377)=v(101)*v(323)+v(100)*v(350)-v(296)*v(98)-v(278)*v(99)
v(376)=v(101)*v(322)+v(100)*v(349)-v(295)*v(98)-v(277)*v(99)
v(375)=-(v(100)*v(249))+v(101)*v(303)-v(330)*v(96)+v(357)*v(99)
v(374)=-(v(100)*v(248))+v(101)*v(302)-v(329)*v(96)+v(356)*v(99)
v(373)=-(v(100)*v(247))+v(101)*v(301)-v(328)*v(96)+v(355)*v(99)
v(372)=-(v(100)*v(246))+v(101)*v(300)-v(327)*v(96)+v(354)*v(99)
v(371)=-(v(100)*v(245))+v(101)*v(299)-v(326)*v(96)+v(353)*v(99)
v(370)=-(v(100)*v(244))+v(101)*v(298)-v(325)*v(96)+v(352)*v(99)
v(369)=-(v(100)*v(243))+v(101)*v(297)-v(324)*v(96)+v(351)*v(99)
v(368)=-(v(100)*v(242))+v(101)*v(296)-v(323)*v(96)+v(350)*v(99)
v(367)=-(v(100)*v(241))+v(101)*v(295)-v(322)*v(96)+v(349)*v(99)
v(366)=-(v(101)*v(267))+v(100)*v(303)-v(357)*v(97)+v(330)*v(99)
v(365)=-(v(101)*v(266))+v(100)*v(302)-v(356)*v(97)+v(329)*v(99)
v(364)=-(v(101)*v(265))+v(100)*v(301)-v(355)*v(97)+v(328)*v(99)
v(363)=-(v(101)*v(264))+v(100)*v(300)-v(354)*v(97)+v(327)*v(99)
v(362)=-(v(101)*v(263))+v(100)*v(299)-v(353)*v(97)+v(326)*v(99)
v(361)=-(v(101)*v(262))+v(100)*v(298)-v(352)*v(97)+v(325)*v(99)
v(360)=-(v(101)*v(261))+v(100)*v(297)-v(351)*v(97)+v(324)*v(99)
v(359)=-(v(101)*v(260))+v(100)*v(296)-v(350)*v(97)+v(323)*v(99)
v(358)=-(v(101)*v(259))+v(100)*v(295)-v(349)*v(97)+v(322)*v(99)
v(132)=-(v(101)*v(97))+v(100)*v(99)
v(127)=-(v(100)*v(96))+v(101)*v(99)
v(122)=v(100)*v(101)-v(98)*v(99)
v(107)=(v(101)*v(101))
v(414)=-(v(103)*v(249))-v(107)*v(267)+v(131)*v(285)+v(303)*v(403)+v(330)*v(405)+v(357)*v(406)-v(339)*v(96)-v(393)*v(97)&
&+v(321)*v(98)
v(413)=-(v(103)*v(248))-v(107)*v(266)+v(131)*v(284)+v(302)*v(403)+v(329)*v(405)+v(356)*v(406)-v(338)*v(96)-v(392)*v(97)&
&+v(320)*v(98)
v(412)=-(v(103)*v(247))-v(107)*v(265)+v(131)*v(283)+v(301)*v(403)+v(328)*v(405)+v(355)*v(406)-v(337)*v(96)-v(391)*v(97)&
&+v(319)*v(98)
v(411)=-(v(103)*v(246))-v(107)*v(264)+v(131)*v(282)+v(300)*v(403)+v(327)*v(405)+v(354)*v(406)-v(336)*v(96)-v(390)*v(97)&
&+v(318)*v(98)
v(410)=-(v(103)*v(245))-v(107)*v(263)+v(131)*v(281)+v(299)*v(403)+v(326)*v(405)+v(353)*v(406)-v(335)*v(96)-v(389)*v(97)&
&+v(317)*v(98)
v(409)=-(v(103)*v(244))-v(107)*v(262)+v(131)*v(280)+v(298)*v(403)+v(325)*v(405)+v(352)*v(406)-v(334)*v(96)-v(388)*v(97)&
&+v(316)*v(98)
v(408)=-(v(103)*v(243))-v(107)*v(261)+v(131)*v(279)+v(297)*v(403)+v(324)*v(405)+v(351)*v(406)-v(333)*v(96)-v(387)*v(97)&
&+v(315)*v(98)
v(407)=-(v(103)*v(242))-v(107)*v(260)+v(131)*v(278)+v(296)*v(403)+v(323)*v(405)+v(350)*v(406)-v(332)*v(96)-v(386)*v(97)&
&+v(314)*v(98)
v(404)=-(v(103)*v(241))-v(107)*v(259)+v(131)*v(277)+v(295)*v(403)+v(322)*v(405)+v(349)*v(406)-v(331)*v(96)-v(385)*v(97)&
&+v(313)*v(98)
v(126)=-v(107)+v(96)*v(98)
v(104)=-(v(103)*v(96))-v(107)*v(97)+v(131)*v(98)+v(403)*v(99)
v(846)=v(122)/v(104)
v(477)=1d0/v(104)**0.23333333333333334d1
v(838)=(-4d0/3d0)*v(477)
v(485)=v(414)*v(838)
v(484)=v(413)*v(838)
v(483)=v(412)*v(838)
v(482)=v(411)*v(838)
v(481)=v(410)*v(838)
v(480)=v(409)*v(838)
v(479)=v(408)*v(838)
v(478)=v(407)*v(838)
v(476)=v(404)*v(838)
v(467)=1d0/v(104)**2
v(475)=-(v(414)*v(467))
v(474)=-(v(413)*v(467))
v(473)=-(v(412)*v(467))
v(472)=-(v(411)*v(467))
v(471)=-(v(410)*v(467))
v(470)=-(v(409)*v(467))
v(469)=-(v(408)*v(467))
v(468)=-(v(407)*v(467))
v(466)=-(v(404)*v(467))
v(426)=1d0/v(104)**0.13333333333333333d1
v(842)=mpar(1)*v(426)
v(839)=-v(426)/3d0
v(435)=v(414)*v(839)
v(434)=v(413)*v(839)
v(433)=v(412)*v(839)
v(432)=v(411)*v(839)
v(431)=v(410)*v(839)
v(430)=v(409)*v(839)
v(429)=v(408)*v(839)
v(428)=v(407)*v(839)
v(427)=v(404)*v(839)
v(415)=sqrt(v(104))
v(840)=mpar(2)*(1d0-1d0/(2d0*v(415)))
v(425)=v(414)*v(840)
v(424)=v(413)*v(840)
v(423)=v(412)*v(840)
v(422)=v(411)*v(840)
v(421)=v(410)*v(840)
v(420)=v(409)*v(840)
v(419)=v(408)*v(840)
v(418)=v(407)*v(840)
v(416)=v(404)*v(840)
v(111)=mpar(2)*(v(104)-v(415))
v(110)=1d0/v(104)**0.3333333333333333d0
v(841)=mpar(1)*v(110)
v(465)=v(425)+mpar(1)*(v(110)*(v(258)+(2d0/3d0)*v(267)+v(294))+v(435)*v(456))
v(464)=v(424)+mpar(1)*(v(110)*(v(257)+(2d0/3d0)*v(266)+v(293))+v(434)*v(456))
v(463)=v(423)+mpar(1)*(v(110)*(v(256)+(2d0/3d0)*v(265)+v(292))+v(433)*v(456))
v(462)=v(422)+mpar(1)*(v(110)*(v(255)+(2d0/3d0)*v(264)+v(291))+v(432)*v(456))
v(461)=v(421)+mpar(1)*(v(110)*(v(254)+(2d0/3d0)*v(263)+v(290))+v(431)*v(456))
v(460)=v(420)+mpar(1)*(v(110)*(v(253)+(2d0/3d0)*v(262)+v(289))+v(430)*v(456))
v(459)=v(419)+mpar(1)*(v(110)*(v(252)+(2d0/3d0)*v(261)+v(288))+v(429)*v(456))
v(458)=v(418)+mpar(1)*(v(110)*(v(251)+(2d0/3d0)*v(260)+v(287))+v(428)*v(456))
v(457)=v(416)+mpar(1)*(v(110)*(v(250)+(2d0/3d0)*v(259)+v(286))+v(427)*v(456))
v(455)=v(425)+mpar(1)*(v(110)*(v(258)+v(276)+(2d0/3d0)*v(285))+v(435)*v(446))
v(454)=v(424)+mpar(1)*(v(110)*(v(257)+v(275)+(2d0/3d0)*v(284))+v(434)*v(446))
v(453)=v(423)+mpar(1)*(v(110)*(v(256)+v(274)+(2d0/3d0)*v(283))+v(433)*v(446))
v(452)=v(422)+mpar(1)*(v(110)*(v(255)+v(273)+(2d0/3d0)*v(282))+v(432)*v(446))
v(451)=v(421)+mpar(1)*(v(110)*(v(254)+v(272)+(2d0/3d0)*v(281))+v(431)*v(446))
v(450)=v(420)+mpar(1)*(v(110)*(v(253)+v(271)+(2d0/3d0)*v(280))+v(430)*v(446))
v(449)=v(419)+mpar(1)*(v(110)*(v(252)+v(270)+(2d0/3d0)*v(279))+v(429)*v(446))
v(448)=v(418)+mpar(1)*(v(110)*(v(251)+v(269)+(2d0/3d0)*v(278))+v(428)*v(446))
v(447)=v(416)+mpar(1)*(v(110)*(v(250)+v(268)+(2d0/3d0)*v(277))+v(427)*v(446))
v(445)=v(425)+mpar(1)*(v(110)*((2d0/3d0)*v(249)+v(276)+v(294))+v(435)*v(436))
v(444)=v(424)+mpar(1)*(v(110)*((2d0/3d0)*v(248)+v(275)+v(293))+v(434)*v(436))
v(443)=v(423)+mpar(1)*(v(110)*((2d0/3d0)*v(247)+v(274)+v(292))+v(433)*v(436))
v(442)=v(422)+mpar(1)*(v(110)*((2d0/3d0)*v(246)+v(273)+v(291))+v(432)*v(436))
v(441)=v(421)+mpar(1)*(v(110)*((2d0/3d0)*v(245)+v(272)+v(290))+v(431)*v(436))
v(440)=v(420)+mpar(1)*(v(110)*((2d0/3d0)*v(244)+v(271)+v(289))+v(430)*v(436))
v(439)=v(419)+mpar(1)*(v(110)*((2d0/3d0)*v(243)+v(270)+v(288))+v(429)*v(436))
v(438)=v(418)+mpar(1)*(v(110)*((2d0/3d0)*v(242)+v(269)+v(287))+v(428)*v(436))
v(437)=v(416)+mpar(1)*(v(110)*((2d0/3d0)*v(241)+v(268)+v(286))+v(427)*v(436))
v(133)=v(111)+v(436)*v(841)
v(848)=v(120)*v(133)
v(845)=v(132)*v(133)
v(128)=v(111)+v(446)*v(841)
v(844)=v(128)*v(131)
v(843)=v(127)*v(128)
v(123)=v(111)+v(456)*v(841)
v(847)=v(123)*v(126)
v(521)=mpar(1)*(v(303)*v(426)+v(485)*v(99))
v(520)=mpar(1)*(v(302)*v(426)+v(484)*v(99))
v(519)=mpar(1)*(v(301)*v(426)+v(483)*v(99))
v(518)=mpar(1)*(v(300)*v(426)+v(482)*v(99))
v(517)=mpar(1)*(v(299)*v(426)+v(481)*v(99))
v(516)=mpar(1)*(v(298)*v(426)+v(480)*v(99))
v(515)=mpar(1)*(v(297)*v(426)+v(479)*v(99))
v(514)=mpar(1)*(v(296)*v(426)+v(478)*v(99))
v(513)=mpar(1)*(v(295)*v(426)+v(476)*v(99))
v(512)=mpar(1)*(v(366)*v(426)+v(132)*v(485))
v(511)=mpar(1)*(v(365)*v(426)+v(132)*v(484))
v(510)=mpar(1)*(v(364)*v(426)+v(132)*v(483))
v(509)=mpar(1)*(v(363)*v(426)+v(132)*v(482))
v(508)=mpar(1)*(v(362)*v(426)+v(132)*v(481))
v(507)=mpar(1)*(v(361)*v(426)+v(132)*v(480))
v(506)=mpar(1)*(v(360)*v(426)+v(132)*v(479))
v(505)=mpar(1)*(v(359)*v(426)+v(132)*v(478))
v(504)=mpar(1)*(v(358)*v(426)+v(132)*v(476))
v(503)=mpar(1)*(v(330)*v(426)+v(100)*v(485))
v(502)=mpar(1)*(v(329)*v(426)+v(100)*v(484))
v(501)=mpar(1)*(v(328)*v(426)+v(100)*v(483))
v(500)=mpar(1)*(v(327)*v(426)+v(100)*v(482))
v(499)=mpar(1)*(v(326)*v(426)+v(100)*v(481))
v(498)=mpar(1)*(v(325)*v(426)+v(100)*v(480))
v(497)=mpar(1)*(v(324)*v(426)+v(100)*v(479))
v(496)=mpar(1)*(v(323)*v(426)+v(100)*v(478))
v(495)=mpar(1)*(v(322)*v(426)+v(100)*v(476))
v(494)=mpar(1)*(v(357)*v(426)+v(101)*v(485))
v(493)=mpar(1)*(v(356)*v(426)+v(101)*v(484))
v(492)=mpar(1)*(v(355)*v(426)+v(101)*v(483))
v(491)=mpar(1)*(v(354)*v(426)+v(101)*v(482))
v(490)=mpar(1)*(v(353)*v(426)+v(101)*v(481))
v(489)=mpar(1)*(v(352)*v(426)+v(101)*v(480))
v(488)=mpar(1)*(v(351)*v(426)+v(101)*v(479))
v(487)=mpar(1)*(v(350)*v(426)+v(101)*v(478))
v(486)=mpar(1)*(v(349)*v(426)+v(101)*v(476))
v(130)=v(101)*v(842)
v(125)=v(100)*v(842)
v(593)=v(130)*v(384)+v(125)*v(402)+(v(128)*v(375)+v(127)*v(455))/v(104)+v(122)*v(494)+v(126)*v(503)+v(475)*v(843)
v(592)=v(130)*v(383)+v(125)*v(401)+(v(128)*v(374)+v(127)*v(454))/v(104)+v(122)*v(493)+v(126)*v(502)+v(474)*v(843)
v(591)=v(130)*v(382)+v(125)*v(400)+(v(128)*v(373)+v(127)*v(453))/v(104)+v(122)*v(492)+v(126)*v(501)+v(473)*v(843)
v(590)=v(130)*v(381)+v(125)*v(399)+(v(128)*v(372)+v(127)*v(452))/v(104)+v(122)*v(491)+v(126)*v(500)+v(472)*v(843)
v(589)=v(130)*v(380)+v(125)*v(398)+(v(128)*v(371)+v(127)*v(451))/v(104)+v(122)*v(490)+v(126)*v(499)+v(471)*v(843)
v(588)=v(130)*v(379)+v(125)*v(397)+(v(128)*v(370)+v(127)*v(450))/v(104)+v(122)*v(489)+v(126)*v(498)+v(470)*v(843)
v(587)=v(130)*v(378)+v(125)*v(396)+(v(128)*v(369)+v(127)*v(449))/v(104)+v(122)*v(488)+v(126)*v(497)+v(469)*v(843)
v(586)=v(130)*v(377)+v(125)*v(395)+(v(128)*v(368)+v(127)*v(448))/v(104)+v(122)*v(487)+v(126)*v(496)+v(468)*v(843)
v(585)=v(130)*v(376)+v(125)*v(394)+(v(128)*v(367)+v(127)*v(447))/v(104)+v(122)*v(486)+v(126)*v(495)+v(466)*v(843)
v(530)=v(125)*v(375)+v(127)*v(503)
v(529)=v(125)*v(374)+v(127)*v(502)
v(528)=v(125)*v(373)+v(127)*v(501)
v(527)=v(125)*v(372)+v(127)*v(500)
v(526)=v(125)*v(371)+v(127)*v(499)
v(525)=v(125)*v(370)+v(127)*v(498)
v(524)=v(125)*v(369)+v(127)*v(497)
v(523)=v(125)*v(368)+v(127)*v(496)
v(522)=v(125)*v(367)+v(127)*v(495)
v(121)=v(132)*v(842)
v(539)=v(121)*v(357)+v(101)*v(512)
v(575)=(v(128)*v(321)+v(131)*v(455))/v(104)+v(530)+v(539)+v(475)*v(844)
v(538)=v(121)*v(356)+v(101)*v(511)
v(574)=(v(128)*v(320)+v(131)*v(454))/v(104)+v(529)+v(538)+v(474)*v(844)
v(537)=v(121)*v(355)+v(101)*v(510)
v(573)=(v(128)*v(319)+v(131)*v(453))/v(104)+v(528)+v(537)+v(473)*v(844)
v(536)=v(121)*v(354)+v(101)*v(509)
v(572)=(v(128)*v(318)+v(131)*v(452))/v(104)+v(527)+v(536)+v(472)*v(844)
v(535)=v(121)*v(353)+v(101)*v(508)
v(571)=(v(128)*v(317)+v(131)*v(451))/v(104)+v(526)+v(535)+v(471)*v(844)
v(534)=v(121)*v(352)+v(101)*v(507)
v(570)=(v(128)*v(316)+v(131)*v(450))/v(104)+v(525)+v(534)+v(470)*v(844)
v(533)=v(121)*v(351)+v(101)*v(506)
v(569)=(v(128)*v(315)+v(131)*v(449))/v(104)+v(524)+v(533)+v(469)*v(844)
v(532)=v(121)*v(350)+v(101)*v(505)
v(568)=(v(128)*v(314)+v(131)*v(448))/v(104)+v(523)+v(532)+v(468)*v(844)
v(531)=v(121)*v(349)+v(101)*v(504)
v(567)=(v(128)*v(313)+v(131)*v(447))/v(104)+v(522)+v(531)+v(466)*v(844)
v(119)=v(842)*v(99)
v(632)=v(130)*v(321)+v(119)*v(375)+(v(133)*v(366)+v(132)*v(445))/v(104)+v(131)*v(494)+v(127)*v(521)+v(475)*v(845)
v(631)=v(130)*v(320)+v(119)*v(374)+(v(133)*v(365)+v(132)*v(444))/v(104)+v(131)*v(493)+v(127)*v(520)+v(474)*v(845)
v(630)=v(130)*v(319)+v(119)*v(373)+(v(133)*v(364)+v(132)*v(443))/v(104)+v(131)*v(492)+v(127)*v(519)+v(473)*v(845)
v(629)=v(130)*v(318)+v(119)*v(372)+(v(133)*v(363)+v(132)*v(442))/v(104)+v(131)*v(491)+v(127)*v(518)+v(472)*v(845)
v(628)=v(130)*v(317)+v(119)*v(371)+(v(133)*v(362)+v(132)*v(441))/v(104)+v(131)*v(490)+v(127)*v(517)+v(471)*v(845)
v(627)=v(130)*v(316)+v(119)*v(370)+(v(133)*v(361)+v(132)*v(440))/v(104)+v(131)*v(489)+v(127)*v(516)+v(470)*v(845)
v(626)=v(130)*v(315)+v(119)*v(369)+(v(133)*v(360)+v(132)*v(439))/v(104)+v(131)*v(488)+v(127)*v(515)+v(469)*v(845)
v(625)=v(130)*v(314)+v(119)*v(368)+(v(133)*v(359)+v(132)*v(438))/v(104)+v(131)*v(487)+v(127)*v(514)+v(468)*v(845)
v(624)=v(130)*v(313)+v(119)*v(367)+(v(133)*v(358)+v(132)*v(437))/v(104)+v(131)*v(486)+v(127)*v(513)+v(466)*v(845)
v(584)=v(121)*v(330)+v(119)*v(348)+v(123)*(v(384)/v(104)+v(122)*v(475))+v(100)*v(512)+v(120)*v(521)+v(465)*v(846)
v(583)=v(121)*v(329)+v(119)*v(347)+v(123)*(v(383)/v(104)+v(122)*v(474))+v(100)*v(511)+v(120)*v(520)+v(464)*v(846)
v(582)=v(121)*v(328)+v(119)*v(346)+v(123)*(v(382)/v(104)+v(122)*v(473))+v(100)*v(510)+v(120)*v(519)+v(463)*v(846)
v(581)=v(121)*v(327)+v(119)*v(345)+v(123)*(v(381)/v(104)+v(122)*v(472))+v(100)*v(509)+v(120)*v(518)+v(462)*v(846)
v(580)=v(121)*v(326)+v(119)*v(344)+v(123)*(v(380)/v(104)+v(122)*v(471))+v(100)*v(508)+v(120)*v(517)+v(461)*v(846)
v(579)=v(121)*v(325)+v(119)*v(343)+v(123)*(v(379)/v(104)+v(122)*v(470))+v(100)*v(507)+v(120)*v(516)+v(460)*v(846)
v(578)=v(121)*v(324)+v(119)*v(342)+v(123)*(v(378)/v(104)+v(122)*v(469))+v(100)*v(506)+v(120)*v(515)+v(459)*v(846)
v(577)=v(121)*v(323)+v(119)*v(341)+v(123)*(v(377)/v(104)+v(122)*v(468))+v(100)*v(505)+v(120)*v(514)+v(458)*v(846)
v(576)=v(121)*v(322)+v(119)*v(340)+v(123)*(v(376)/v(104)+v(122)*v(466))+v(100)*v(504)+v(120)*v(513)+v(457)*v(846)
v(548)=v(119)*v(384)+v(122)*v(521)
v(566)=(v(123)*v(402)+v(126)*v(465))/v(104)+v(530)+v(548)+v(475)*v(847)
v(557)=(v(133)*v(348)+v(120)*v(445))/v(104)+v(539)+v(548)+v(475)*v(848)
v(547)=v(119)*v(383)+v(122)*v(520)
v(565)=(v(123)*v(401)+v(126)*v(464))/v(104)+v(529)+v(547)+v(474)*v(847)
v(556)=(v(133)*v(347)+v(120)*v(444))/v(104)+v(538)+v(547)+v(474)*v(848)
v(546)=v(119)*v(382)+v(122)*v(519)
v(564)=(v(123)*v(400)+v(126)*v(463))/v(104)+v(528)+v(546)+v(473)*v(847)
v(555)=(v(133)*v(346)+v(120)*v(443))/v(104)+v(537)+v(546)+v(473)*v(848)
v(545)=v(119)*v(381)+v(122)*v(518)
v(563)=(v(123)*v(399)+v(126)*v(462))/v(104)+v(527)+v(545)+v(472)*v(847)
v(554)=(v(133)*v(345)+v(120)*v(442))/v(104)+v(536)+v(545)+v(472)*v(848)
v(544)=v(119)*v(380)+v(122)*v(517)
v(562)=(v(123)*v(398)+v(126)*v(461))/v(104)+v(526)+v(544)+v(471)*v(847)
v(553)=(v(133)*v(344)+v(120)*v(441))/v(104)+v(535)+v(544)+v(471)*v(848)
v(543)=v(119)*v(379)+v(122)*v(516)
v(561)=(v(123)*v(397)+v(126)*v(460))/v(104)+v(525)+v(543)+v(470)*v(847)
v(552)=(v(133)*v(343)+v(120)*v(440))/v(104)+v(534)+v(543)+v(470)*v(848)
v(542)=v(119)*v(378)+v(122)*v(515)
v(560)=(v(123)*v(396)+v(126)*v(459))/v(104)+v(524)+v(542)+v(469)*v(847)
v(551)=(v(133)*v(342)+v(120)*v(439))/v(104)+v(533)+v(542)+v(469)*v(848)
v(541)=v(119)*v(377)+v(122)*v(514)
v(559)=(v(123)*v(395)+v(126)*v(458))/v(104)+v(523)+v(541)+v(468)*v(847)
v(550)=(v(133)*v(341)+v(120)*v(438))/v(104)+v(532)+v(541)+v(468)*v(848)
v(540)=v(119)*v(376)+v(122)*v(513)
v(558)=(v(123)*v(394)+v(126)*v(457))/v(104)+v(522)+v(540)+v(466)*v(847)
v(549)=(v(133)*v(340)+v(120)*v(437))/v(104)+v(531)+v(540)+v(466)*v(848)
v(114)=v(125)*v(127)
v(113)=v(101)*v(121)
v(106)=v(119)*v(122)
v(105)=v(106)+v(113)+v(848)/v(104)
v(112)=v(106)+v(114)+v(847)/v(104)
v(118)=v(113)+v(114)+v(844)/v(104)
v(124)=v(119)*v(120)+v(100)*v(121)+v(123)*v(846)
v(129)=v(125)*v(126)+v(122)*v(130)+v(843)/v(104)
v(613)=v(124)*v(233)+v(112)*v(235)+v(129)*v(240)
v(616)=v(613)+v(559)*v(79)+v(586)*v(88)+v(577)*v(94)
v(609)=v(124)*v(232)+v(112)*v(237)+v(129)*v(239)
v(622)=v(609)+v(565)*v(79)+v(592)*v(88)+v(583)*v(94)
v(605)=v(124)*v(234)+v(112)*v(236)+v(129)*v(238)
v(619)=v(605)+v(562)*v(79)+v(589)*v(88)+v(580)*v(94)
v(600)=v(605)+v(582)*v(78)+v(564)*v(84)+v(591)*v(93)
v(597)=v(613)+v(579)*v(78)+v(561)*v(84)+v(588)*v(93)
v(594)=v(609)+v(576)*v(78)+v(558)*v(84)+v(585)*v(93)
v(149)=v(124)*v(78)+v(112)*v(84)+v(129)*v(93)
v(145)=v(129)*v(80)+v(124)*v(92)+v(112)*v(95)
v(141)=v(112)*v(79)+v(129)*v(88)+v(124)*v(94)
v(134)=v(119)*v(127)+v(130)*v(131)+v(845)/v(104)
v(673)=v(105)*v(233)+v(124)*v(235)+v(134)*v(240)
v(685)=v(673)+v(577)*v(79)+v(625)*v(88)+v(550)*v(94)
v(669)=v(105)*v(232)+v(124)*v(237)+v(134)*v(239)
v(691)=v(669)+v(583)*v(79)+v(631)*v(88)+v(556)*v(94)
v(665)=v(105)*v(234)+v(124)*v(236)+v(134)*v(238)
v(688)=v(665)+v(580)*v(79)+v(628)*v(88)+v(553)*v(94)
v(661)=v(134)*v(233)+v(129)*v(235)+v(118)*v(240)
v(676)=v(661)+v(586)*v(79)+v(568)*v(88)+v(625)*v(94)
v(657)=v(134)*v(232)+v(129)*v(237)+v(118)*v(239)
v(682)=v(657)+v(592)*v(79)+v(574)*v(88)+v(631)*v(94)
v(653)=v(134)*v(234)+v(129)*v(236)+v(118)*v(238)
v(679)=v(653)+v(589)*v(79)+v(571)*v(88)+v(628)*v(94)
v(648)=v(653)+v(630)*v(78)+v(591)*v(84)+v(573)*v(93)
v(645)=v(661)+v(627)*v(78)+v(588)*v(84)+v(570)*v(93)
v(642)=v(657)+v(624)*v(78)+v(585)*v(84)+v(567)*v(93)
v(639)=v(665)+v(555)*v(78)+v(582)*v(84)+v(630)*v(93)
v(636)=v(673)+v(552)*v(78)+v(579)*v(84)+v(627)*v(93)
v(633)=v(669)+v(549)*v(78)+v(576)*v(84)+v(624)*v(93)
v(151)=v(105)*v(78)+v(124)*v(84)+v(134)*v(93)
v(150)=v(134)*v(78)+v(129)*v(84)+v(118)*v(93)
v(147)=v(118)*v(80)+v(134)*v(92)+v(129)*v(95)
v(146)=v(134)*v(80)+v(105)*v(92)+v(124)*v(95)
v(143)=v(129)*v(79)+v(118)*v(88)+v(134)*v(94)
v(142)=v(124)*v(79)+v(134)*v(88)+v(105)*v(94)
v(741)=v(142)*v(232)+v(141)*v(237)+v(143)*v(239)
v(736)=v(142)*v(234)+v(141)*v(236)+v(143)*v(238)
v(714)=v(142)*v(233)+v(141)*v(235)+v(143)*v(240)
v(144)=v(142)*v(78)+v(141)*v(84)+v(143)*v(93)
v(864)=v(144)*v(788)
v(148)=v(145)*v(79)+v(147)*v(88)+v(146)*v(94)
v(866)=v(148)*v(788)
v(152)=v(150)*v(80)+v(151)*v(92)+v(149)*v(95)
v(865)=v(152)*v(788)
v(160)=1d0/(statev(16)*v(179)+statev(19)*v(180)+v(178)*v(796))**2
v(167)=-(v(160)*((v(178)*v(178))+(v(179)*v(179))+(v(180)*v(180))))
v(165)=1d0/v(160)**0.3333333333333333d0
v(851)=-(mpar(13)*v(165))
v(852)=v(160)*v(851)
v(164)=v(160)*((v(172)*v(172))+(v(173)*v(173))+(v(176)*v(176)))
v(169)=-v(164)/3d0
v(163)=-(v(160)*((v(171)*v(171))+(v(174)*v(174))+(v(175)*v(175))))
v(168)=v(163)/3d0
v(162)=v(167)/3d0
v(183)=1d0/(statev(25)*v(202)+statev(28)*v(203)+v(201)*v(808))**2
v(190)=-(v(183)*((v(201)*v(201))+(v(202)*v(202))+(v(203)*v(203))))
v(188)=1d0/v(183)**0.3333333333333333d0
v(849)=mpar(15)*v(188)
v(850)=v(183)*v(849)
v(187)=v(183)*((v(195)*v(195))+(v(196)*v(196))+(v(199)*v(199)))
v(192)=-v(187)/3d0
v(186)=-(v(183)*((v(194)*v(194))+(v(197)*v(197))+(v(198)*v(198))))
v(191)=v(186)/3d0
v(185)=v(190)/3d0
v(184)=(v(185)+(2d0/3d0)*v(187)+v(191))*v(849)
v(189)=(v(185)+(-2d0/3d0)*v(186)+v(192))*v(849)
v(200)=(v(194)*v(195)+v(196)*v(197)+v(198)*v(199))*v(850)
v(204)=(v(194)*v(201)+v(198)*v(202)+v(197)*v(203))*v(850)
v(205)=(v(195)*v(201)+v(199)*v(202)+v(196)*v(203))*v(850)
v(206)=v(133)-v(184)+(v(162)+(2d0/3d0)*v(164)+v(168))*v(851)
v(207)=v(123)-v(189)+(v(162)+(-2d0/3d0)*v(163)+v(169))*v(851)
v(215)=-v(207)/3d0
v(208)=v(128)+((2d0/3d0)*v(190)-v(191)-v(192))*v(849)+((-2d0/3d0)*v(167)+v(168)+v(169))*v(851)
v(214)=-v(208)/3d0
v(209)=v(101)*v(129)-v(200)+(v(171)*v(172)+v(173)*v(174)+v(175)*v(176))*v(852)+v(124)*v(96)+v(112)*v(99)
v(210)=v(100)*v(118)-v(204)+(v(171)*v(178)+v(175)*v(179)+v(174)*v(180))*v(852)+v(129)*v(97)+v(134)*v(99)
v(211)=v(101)*v(105)+v(100)*v(124)-v(205)+(v(172)*v(178)+v(176)*v(179)+v(173)*v(180))*v(852)+v(134)*v(98)
v(212)=(2d0/3d0)*v(206)+v(214)+v(215)
v(213)=-v(206)/3d0
v(220)=(2d0/3d0)*v(208)+v(213)+v(215)
v(218)=(2d0/3d0)*v(207)+v(213)+v(214)
v(225)=0.10000000000000012d1*sqrt(0.1d-14+v(825)**2+v(853)+v(854)+v(855)+v(856)+v(857))
v(228)=(v(858)*v(858))/(statev(34)*v(225)+v(858))**2
v(699)=v(550)*v(859)+v(559)*v(860)+v(568)*v(861)+v(586)*v(862)+v(826)*(v(577)*v(84)+v(625)*v(93))
v(705)=v(551)*v(859)+v(560)*v(860)+v(569)*v(861)+v(587)*v(862)+v(826)*(v(578)*v(84)+v(626)*v(93))
v(718)=v(553)*v(859)+v(562)*v(860)+v(571)*v(861)+v(589)*v(862)+v(826)*(v(580)*v(84)+v(628)*v(93))
v(725)=v(554)*v(859)+v(563)*v(860)+v(572)*v(861)+v(590)*v(862)+v(826)*(v(581)*v(84)+v(629)*v(93))
v(781)=v(863)*(Fnew(8)*(v(741)+v(78)*(v(576)*v(79)+v(624)*v(88)+v(549)*v(94))+v(84)*(v(558)*v(79)+v(585)*v(88)+v(576)*v&
&(94))+v(93)*(v(585)*v(79)+v(567)*v(88)+v(624)*v(94)))+Fnew(2)*(v(714)+v(78)*(v(579)*v(79)+v(627)*v(88)+v(552)*v(94))+v&
&(84)*(v(561)*v(79)+v(588)*v(88)+v(579)*v(94))+v(93)*(v(588)*v(79)+v(570)*v(88)+v(627)*v(94)))+Fnew(5)*(v(736)+v(78)*(v&
&(582)*v(79)+v(630)*v(88)+v(555)*v(94))+v(84)*(v(564)*v(79)+v(591)*v(88)+v(582)*v(94))+v(93)*(v(591)*v(79)+v(573)*v(88)&
&+v(630)*v(94))))
v(780)=v(863)*(Fnew(6)*(v(642)*v(80)+v(633)*v(92)+v(594)*v(95))+Fnew(9)*(v(645)*v(80)+v(636)*v(92)+v(597)*v(95))+Fnew(3&
&)*(v(648)*v(80)+v(639)*v(92)+v(600)*v(95)))
v(740)=v(556)*v(859)+v(565)*v(860)+v(574)*v(861)+v(592)*v(862)+v(826)*(v(583)*v(84)+v(631)*v(93))
v(776)=v(863)*(Fnew(4)*(v(685)*v(78)+v(616)*v(84)+v(676)*v(93))+Fnew(7)*(v(688)*v(78)+v(619)*v(84)+v(679)*v(93))+Fnew(1&
&)*(v(691)*v(78)+v(622)*v(84)+v(682)*v(93)))
v(748)=v(557)*v(859)+v(566)*v(860)+v(575)*v(861)+v(593)*v(862)+v(826)*(v(584)*v(84)+v(632)*v(93))
v(759)=(Fnew(4)*v(699)+Fnew(7)*v(718)+Fnew(1)*v(740))*v(863)
v(760)=-v(759)/4d0
v(761)=(Fnew(7)*v(705)+Fnew(1)*v(725)+Fnew(4)*v(748))*v(863)
v(762)=-v(761)/4d0
v(763)=(Fnew(9)*v(699)+Fnew(3)*v(718)+Fnew(6)*v(740))*v(863)
v(768)=-v(763)/4d0
sigma(1)=v(788)*(v(151)*v(78)+v(149)*v(84)+v(150)*v(93))
sigma(2)=v(788)*(v(141)*v(79)+v(143)*v(88)+v(142)*v(94))
sigma(3)=v(788)*(v(147)*v(80)+v(146)*v(92)+v(145)*v(95))
sigma(4)=v(864)
sigma(5)=v(865)
sigma(6)=v(866)
ddsdde(1,1)=v(788)*(Fnew(1)*(v(151)*v(232)+v(149)*v(237)+v(150)*v(239)+v(633)*v(78)+v(594)*v(84)+v(642)*v(93))+Fnew(4)*&
&(v(151)*v(233)+v(149)*v(235)+v(150)*v(240)+v(636)*v(78)+v(597)*v(84)+v(645)*v(93))+Fnew(7)*(v(151)*v(234)+v(149)*v(236)&
&+v(150)*v(238)+v(639)*v(78)+v(600)*v(84)+v(648)*v(93)))
ddsdde(1,2)=(Fnew(2)*v(699)+Fnew(5)*v(718)+Fnew(8)*v(740))*v(788)
ddsdde(1,3)=(Fnew(3)*v(705)+Fnew(6)*v(725)+Fnew(9)*v(748))*v(788)
ddsdde(1,4)=v(760)
ddsdde(1,5)=v(762)
ddsdde(1,6)=v(763)/2d0
ddsdde(2,2)=v(788)*(Fnew(2)*(v(714)+v(616)*v(79)+v(676)*v(88)+v(685)*v(94))+Fnew(5)*(v(736)+v(619)*v(79)+v(679)*v(88)+v&
&(688)*v(94))+Fnew(8)*(v(741)+v(622)*v(79)+v(682)*v(88)+v(691)*v(94)))
ddsdde(2,3)=v(788)*(Fnew(3)*(v(94)*(v(578)*v(79)+v(626)*v(88)+v(551)*v(94))+v(79)*(v(560)*v(79)+v(587)*v(88)+v(578)*v&
&(94))+v(88)*(v(587)*v(79)+v(569)*v(88)+v(626)*v(94)))+Fnew(6)*(v(94)*(v(581)*v(79)+v(629)*v(88)+v(554)*v(94))+v(79)*(v&
&(563)*v(79)+v(590)*v(88)+v(581)*v(94))+v(88)*(v(590)*v(79)+v(572)*v(88)+v(629)*v(94)))+Fnew(9)*(v(94)*(v(584)*v(79)+v&
&(632)*v(88)+v(557)*v(94))+v(79)*(v(566)*v(79)+v(593)*v(88)+v(584)*v(94))+v(88)*(v(593)*v(79)+v(575)*v(88)+v(632)*v(94))&
&))
ddsdde(2,4)=v(760)
ddsdde(2,5)=v(761)/2d0
ddsdde(2,6)=v(768)
ddsdde(3,3)=v(788)*(Fnew(3)*(v(146)*v(234)+v(145)*v(236)+v(147)*v(238)+v(95)*(v(605)+v(587)*v(80)+v(578)*v(92)+v(560)*v&
&(95))+v(92)*(v(665)+v(626)*v(80)+v(551)*v(92)+v(578)*v(95))+v(80)*(v(653)+v(569)*v(80)+v(626)*v(92)+v(587)*v(95)))+Fnew&
&(6)*(v(146)*v(232)+v(145)*v(237)+v(147)*v(239)+v(95)*(v(609)+v(590)*v(80)+v(581)*v(92)+v(563)*v(95))+v(92)*(v(669)+v&
&(629)*v(80)+v(554)*v(92)+v(581)*v(95))+v(80)*(v(657)+v(572)*v(80)+v(629)*v(92)+v(590)*v(95)))+Fnew(9)*(v(146)*v(233)+v&
&(145)*v(235)+v(147)*v(240)+v(95)*(v(613)+v(593)*v(80)+v(584)*v(92)+v(566)*v(95))+v(92)*(v(673)+v(632)*v(80)+v(557)*v(92&
&)+v(584)*v(95))+v(80)*(v(661)+v(575)*v(80)+v(632)*v(92)+v(593)*v(95))))
ddsdde(3,4)=v(759)/2d0
ddsdde(3,5)=v(762)
ddsdde(3,6)=v(768)
ddsdde(4,4)=(v(776)+v(781))/4d0
ddsdde(4,5)=v(866)/2d0
ddsdde(4,6)=v(865)/2d0
ddsdde(5,5)=(v(776)+v(780))/4d0
ddsdde(5,6)=v(864)/2d0
ddsdde(6,6)=(v(780)+v(781))/4d0
DO i01=2,6
  DO i02=1,i01-1
    ddsdde(i01,i02)=ddsdde(i02,i01)
  ENDDO
ENDDO
yielding=-statev(36)-v(858)+sqrt(0.15d1*((2d0*v(209)**2+2d0*v(210)**2+2d0*v(211)**2+v(212)**2+v(218)**2+v(220)**2)*v&
&(228)+((1d0-v(228))*(2d0*statev(31)*v(209)+2d0*statev(33)*v(210)+2d0*statev(32)*v(211)+statev(29)*v(212)+statev(30)*v&
&(218)+v(220)*v(825))**2)/v(225)**2))
xguess(1)=statev(10)
xguess(2)=v(212)
xguess(3)=v(218)
xguess(4)=v(209)
xguess(5)=v(211)
xguess(6)=v(210)
xguess(7)=v(184)
xguess(8)=v(189)
xguess(9)=v(200)
xguess(10)=v(205)
xguess(11)=v(204)
xguess(12)=statev(29)
xguess(13)=statev(30)
xguess(14)=statev(31)
xguess(15)=statev(32)
xguess(16)=statev(33)
END