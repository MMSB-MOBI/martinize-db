"use strict";
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var __generator = (this && this.__generator) || function (thisArg, body) {
    var _ = { label: 0, sent: function() { if (t[0] & 1) throw t[1]; return t[1]; }, trys: [], ops: [] }, f, y, t, g;
    return g = { next: verb(0), "throw": verb(1), "return": verb(2) }, typeof Symbol === "function" && (g[Symbol.iterator] = function() { return this; }), g;
    function verb(n) { return function (v) { return step([n, v]); }; }
    function step(op) {
        if (f) throw new TypeError("Generator is already executing.");
        while (_) try {
            if (f = 1, y && (t = op[0] & 2 ? y["return"] : op[0] ? y["throw"] || ((t = y["return"]) && t.call(y), 0) : y.next) && !(t = t.call(y, op[1])).done) return t;
            if (y = 0, t) op = [op[0] & 2, t.value];
            switch (op[0]) {
                case 0: case 1: t = op; break;
                case 4: _.label++; return { value: op[1], done: false };
                case 5: _.label++; y = op[1]; op = [0]; continue;
                case 7: op = _.ops.pop(); _.trys.pop(); continue;
                default:
                    if (!(t = _.trys, t = t.length > 0 && t[t.length - 1]) && (op[0] === 6 || op[0] === 2)) { _ = 0; continue; }
                    if (op[0] === 3 && (!t || (op[1] > t[0] && op[1] < t[3]))) { _.label = op[1]; break; }
                    if (op[0] === 6 && _.label < t[1]) { _.label = t[1]; t = op; break; }
                    if (t && _.label < t[2]) { _.label = t[2]; _.ops.push(op); break; }
                    if (t[2]) _.ops.pop();
                    _.trys.pop(); continue;
            }
            op = body.call(thisArg, _);
        } catch (e) { op = [6, e]; y = 0; } finally { f = t = 0; }
        if (op[0] & 5) throw op[1]; return { value: op[0] ? op[1] : void 0, done: true };
    }
};
exports.__esModule = true;
var nodemailer_1 = require("nodemailer");
var constants_1 = require("../constants");
var twig_1 = require("twig");
var logger_1 = require("../logger");
exports["default"] = new /** @class */ (function () {
    function Mailer() {
        this.transporter = nodemailer_1["default"].createTransport(constants_1.MAILER_TRANSPORT_SETTINGS);
        this.default_sender = { name: constants_1.DEFAULT_MAILER_NAME, address: constants_1.DEFAULT_MAILER_ADDRESS };
    }
    Mailer.prototype.send = function (send_options, template_name, options) {
        return __awaiter(this, void 0, void 0, function () {
            var file, content;
            return __generator(this, function (_a) {
                switch (_a.label) {
                    case 0:
                        if (!send_options.to) {
                            throw new Error("You must define a mail recipient.");
                        }
                        if (!options.site_url) {
                            options.site_url = constants_1.URLS.SERVER;
                        }
                        if (!options.static_site_url) {
                            options.static_site_url = constants_1.URLS.SERVER;
                        }
                        file = constants_1.TEMPLATE_DIR + template_name + (template_name.endsWith('.twig') ? "" : ".twig");
                        return [4 /*yield*/, new Promise(function (resolve, reject) {
                                // @ts-ignore Incorrect typedef for options
                                twig_1["default"].renderFile(file, options, function (err, res) {
                                    if (err) {
                                        reject(err);
                                    }
                                    resolve(res);
                                });
                            })];
                    case 1:
                        content = _a.sent();
                        send_options.html = content;
                        if (!send_options.from) {
                            send_options.from = this.default_sender;
                            send_options.sender = this.default_sender;
                        }
                        if (!send_options.subject && options.title) {
                            send_options.subject = options.title;
                        }
                        if (constants_1.MAILER_ENFORCE_RECIPIENT) {
                            send_options.to = constants_1.MAILER_ENFORCE_RECIPIENT;
                        }
                        return [2 /*return*/, this.mail(send_options)];
                }
            });
        });
    };
    Mailer.prototype.mail = function (options) {
        return __awaiter(this, void 0, void 0, function () {
            var info, e_1;
            return __generator(this, function (_a) {
                switch (_a.label) {
                    case 0:
                        _a.trys.push([0, 2, , 3]);
                        return [4 /*yield*/, this.transporter.sendMail(options)];
                    case 1:
                        info = _a.sent();
                        logger_1["default"].debug('Sended email:' + info.messageId);
                        return [2 /*return*/, info];
                    case 2:
                        e_1 = _a.sent();
                        logger_1["default"].error('Unable to send email to ' + options.sender + ' / ' + e_1);
                        logger_1["default"].debug(options);
                        logger_1["default"].debug(e_1);
                        throw new Error('Unable to send email.');
                    case 3: return [2 /*return*/];
                }
            });
        });
    };
    return Mailer;
}())();
