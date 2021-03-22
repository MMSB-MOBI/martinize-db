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
exports.SearchWorker = void 0;
var worker_thread_manager_1 = require("worker-thread-manager");
var CouchHelper_1 = require("../Entities/CouchHelper");
var constants_1 = require("../constants");
exports.SearchWorker = new /** @class */ (function () {
    function SearchWorker() {
    }
    Object.defineProperty(SearchWorker.prototype, "pool", {
        get: function () {
            if (!this._pool) {
                this._pool = worker_thread_manager_1["default"].spawn(__dirname + '/worker.js', {
                    spawnerThreshold: constants_1.MAX_REQUEST_PER_WORKER_THRESHOLD,
                    poolLength: constants_1.MAX_POOL_SIZE,
                    // 2 minutes life max without task
                    stopOnNoTask: constants_1.MINUTES_BEFORE_WORKER_KILL * 60 * 1000,
                    workerData: {
                        molecule_collection: CouchHelper_1["default"].MOLECULE_COLLECTION,
                        couch_url: constants_1.URLS.COUCH
                    }
                });
            }
            return this._pool;
        },
        enumerable: false,
        configurable: true
    });
    SearchWorker.prototype.query = function (query, send_all) {
        if (send_all === void 0) { send_all = false; }
        return __awaiter(this, void 0, void 0, function () {
            var task;
            return __generator(this, function (_a) {
                task = {
                    type: 'new_search',
                    query: query,
                    as_all: send_all
                };
                return [2 /*return*/, this.pool.run(task)];
            });
        });
    };
    SearchWorker.prototype.clearCache = function () {
        return this.pool.send({ type: 'clean' });
    };
    return SearchWorker;
}());
exports["default"] = exports.SearchWorker;
