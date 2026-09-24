#pragma once
#include <cstdint>
#include <string>

namespace fbs::caves::detail {
inline bool valid_utf8(const std::string& text) {
    for (std::size_t i=0;i<text.size();) {
        const auto lead=static_cast<unsigned char>(text[i++]);
        if (lead<0x80) continue;
        std::uint32_t value=0,minimum=0;
        unsigned remaining=0;
        if (lead>=0xc2 && lead<=0xdf) {value=lead&0x1f;minimum=0x80;remaining=1;}
        else if (lead>=0xe0 && lead<=0xef) {value=lead&0x0f;minimum=0x800;remaining=2;}
        else if (lead>=0xf0 && lead<=0xf4) {value=lead&0x07;minimum=0x10000;remaining=3;}
        else return false;
        if (text.size()-i<remaining) return false;
        while (remaining--) {
            const auto byte=static_cast<unsigned char>(text[i++]);
            if ((byte&0xc0)!=0x80) return false;
            value=(value<<6)|(byte&0x3f);
        }
        if (value<minimum || value>0x10ffff || (value>=0xd800 && value<=0xdfff)) return false;
    }
    return true;
}
}
